import gzip
import hashlib
import io

import pytest

from normalize_tes import phi_sfs
from normalize_tes import vcf_eligibility
from normalize_tes.vcf_io import (
    COMPRESSED_SUFFIXES,
    _decode_genotype,
    _HashingStream,
    _open_vcf,
    decode_inbred_genotype,
    decode_inbred_genotype_short,
    drain,
    open_vcf,
)


def _read_all_and_drain(handle, buffered):
    text = handle.read()
    drain(buffered)
    return text


def test_open_vcf_hash_matches_raw_bytes_for_plain_text(tmp_path):
    path = tmp_path / "sites.vcf"
    raw = b"##fileformat=VCFv4.2\nchr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0\n"
    path.write_bytes(raw)
    handle, hashing, buffered = open_vcf(path)
    try:
        text = _read_all_and_drain(handle, buffered)
    finally:
        handle.close()
    assert text == raw.decode("utf-8")
    assert hashing.digest.hexdigest() == hashlib.sha256(raw).hexdigest()


@pytest.mark.parametrize("suffix", [".gz", ".bgz", ".bgzf"])
def test_open_vcf_hash_matches_raw_bytes_for_gzip_suffixes(tmp_path, suffix):
    path = tmp_path / f"sites.vcf{suffix}"
    payload = b"##fileformat=VCFv4.2\nchr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0\n"
    compressed = gzip.compress(payload)
    path.write_bytes(compressed)
    handle, hashing, buffered = open_vcf(path)
    try:
        text = _read_all_and_drain(handle, buffered)
    finally:
        handle.close()
    assert text == payload.decode("utf-8")
    # The digest covers the compressed bytes on disk, not the decompressed
    # text, since it authenticates the input file's identity.
    assert hashing.digest.hexdigest() == hashlib.sha256(compressed).hexdigest()


def test_open_vcf_treats_non_gz_suffix_as_uncompressed(tmp_path):
    path = tmp_path / "sites.vcf.txt"
    raw = b"plain text, not gzip\n"
    path.write_bytes(raw)
    handle, hashing, buffered = open_vcf(path)
    try:
        text = _read_all_and_drain(handle, buffered)
    finally:
        handle.close()
    assert text == raw.decode("utf-8")
    assert hashing.digest.hexdigest() == hashlib.sha256(raw).hexdigest()


def test_open_vcf_multi_member_gzip_digest_covers_every_member(tmp_path):
    """A digest that stopped at the first gzip member would miss real bytes.

    ``drain`` exists because the text/gzip layers can stop consuming before
    physical EOF; a multi-member file is the concrete case where that
    matters, since GzipFile transparently concatenates members but the
    hashing stream sees each member's trailing bytes only if the buffered
    reader is drained after the caller is done.
    """
    first = b"chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0\n"
    second = b"chr1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t1\n"
    compressed = gzip.compress(first) + gzip.compress(second)
    path = tmp_path / "sites.vcf.gz"
    path.write_bytes(compressed)
    handle, hashing, buffered = open_vcf(path)
    try:
        text = _read_all_and_drain(handle, buffered)
    finally:
        handle.close()
    assert text == (first + second).decode("utf-8")
    assert hashing.digest.hexdigest() == hashlib.sha256(compressed).hexdigest()
    # Sanity check: the digest is NOT just the first member's compressed
    # bytes, i.e. draining genuinely picked up the second member.
    assert hashing.digest.hexdigest() != hashlib.sha256(gzip.compress(first)).hexdigest()


def test_compressed_suffixes_and_aliases_are_exported():
    assert COMPRESSED_SUFFIXES == (".gz", ".bgz", ".bgzf")
    assert _open_vcf is open_vcf
    assert _HashingStream.__name__ == "_HashingStream"


# --- genotype decoder --------------------------------------------------

# (gt, heterozygous_policy, expected canonical tag, expected phi_sfs short tag)
_DECODE_CASES = [
    ("0", "error", 0, 0),
    ("1", "error", 1, 1),
    ("0|0", "error", 0, 0),
    ("1|1", "error", 1, 1),
    ("0/0", "error", 0, 0),
    ("1/1", "error", 1, 1),
    ("0/0/0", "error", 0, 0),  # triploid inbred homozygous call
    (".", "error", None, None),
    ("./.", "error", None, None),
    ("0/.", "error", None, None),
    (".|0", "error", None, None),
    ("0/1", "error", "heterozygous_genotype", "heterozygous"),
    ("1|0", "error", "heterozygous_genotype", "heterozygous"),
    ("0/1", "missing", None, None),
    ("1|0", "missing", None, None),
    ("2", "error", "non_biallelic_genotype", "non-biallelic"),
    ("0/2", "error", "non_biallelic_genotype", "non-biallelic"),
    ("2/2", "missing", "non_biallelic_genotype", "non-biallelic"),
    ("x", "error", "invalid_genotype", "invalid"),
    ("-1", "error", "invalid_genotype", "invalid"),
    ("0/x", "missing", "invalid_genotype", "invalid"),
]


@pytest.mark.parametrize("gt, heterozygous, canonical, short", _DECODE_CASES)
def test_decode_inbred_genotype_canonical_tags(gt, heterozygous, canonical, short):
    assert decode_inbred_genotype(gt, heterozygous) == canonical


@pytest.mark.parametrize("gt, heterozygous, canonical, short", _DECODE_CASES)
def test_decode_inbred_genotype_short_tags_match_phi_sfs_table(gt, heterozygous, canonical, short):
    result = decode_inbred_genotype_short(gt, heterozygous)
    assert result == short
    # phi_sfs looks up a string result in its own _GENOTYPE_ERRORS table; any
    # tag this function returns must be one of that table's keys so a
    # drop-in caller does not KeyError.
    if isinstance(result, str):
        assert result in phi_sfs._GENOTYPE_ERRORS
    assert _decode_genotype is decode_inbred_genotype_short


@pytest.mark.parametrize("gt, heterozygous, canonical, short", _DECODE_CASES)
def test_decoder_agrees_with_phi_sfs_current_decoder(gt, heterozygous, canonical, short):
    """The equal-M invariant depends on both readers making the identical
    callability decision for every GT string. This compares vcf_io's short-
    tag wrapper against phi_sfs's own (still separately defined) decoder,
    read-only, without importing phi_sfs's private helper as a replacement.
    """
    assert decode_inbred_genotype_short(gt, heterozygous) == phi_sfs._decode_genotype(
        gt, heterozygous
    )


@pytest.mark.parametrize("gt, heterozygous, canonical, short", _DECODE_CASES)
def test_decoder_agrees_with_vcf_eligibility_current_decoder(gt, heterozygous, canonical, short):
    """vcf_eligibility now imports vcf_io.decode_inbred_genotype directly, so
    this is mostly a guard against a future re-divergence if that import is
    ever replaced by a local copy again.
    """
    assert decode_inbred_genotype(gt, heterozygous) == vcf_eligibility.decode_inbred_genotype(
        gt, heterozygous
    )
