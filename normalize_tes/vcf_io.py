#!/usr/bin/env python3
"""Shared low-level VCF reading primitives.

`vcf_eligibility.scan_vcf` and `phi_sfs.read_site_counts` both open the
analysis VCF once, digest its raw bytes while a text/gzip layer decodes it,
and decode each inbred individual's GT string under the same missing/
heterozygous rules. Keeping one copy of that plumbing here means the two
readers cannot silently drift apart: the equal-M invariant between an
eligibility mask and a `phi_sfs` run depends on both making identical
callability decisions for the same GT string.

The two callers report genotype-decode failures under different reason
tags (`vcf_eligibility` uses snake_case report keys, `phi_sfs` uses short
hyphenated tags matched against its own error-message table). Rather than
force one caller to change its public strings, `decode_inbred_genotype`
returns the canonical (snake_case) tag and `decode_inbred_genotype_short`
is a thin wrapper that remaps it to `phi_sfs`'s tags. `_decode_genotype` is
an alias for the short-tag wrapper so `phi_sfs` can import it under its
current private name with no call-site changes.
"""

from __future__ import annotations

import gzip
import hashlib
import io
from pathlib import Path

COMPRESSED_SUFFIXES = (".gz", ".bgz", ".bgzf")


class _HashingStream(io.RawIOBase):
    """Raw byte stream that digests everything read through it."""

    def __init__(self, handle):
        self._handle = handle
        self.digest = hashlib.sha256()

    def readable(self) -> bool:
        return True

    def readinto(self, buffer) -> int:
        count = self._handle.readinto(buffer)
        if count:
            self.digest.update(memoryview(buffer)[:count])
        return count

    def close(self) -> None:
        try:
            self._handle.close()
        finally:
            super().close()


def open_vcf(path: Path):
    """Open a VCF as text over a hashing stream, so one pass yields both.

    Returns `(text_handle, hashing, buffered)`. The caller reads records
    from `text_handle` and, once exhausted, should call `drain(buffered)`
    before reading `hashing.digest`, so the digest always covers every
    physical input byte -- including bytes a text decoder or gzip member
    boundary left unread (see `drain`).
    """
    path = Path(path)
    hashing = _HashingStream(path.open("rb"))
    buffered = io.BufferedReader(hashing, buffer_size=1 << 20)
    stream = (
        gzip.GzipFile(fileobj=buffered)
        if path.suffix.lower() in COMPRESSED_SUFFIXES
        else buffered
    )
    return io.TextIOWrapper(stream, encoding="utf-8"), hashing, buffered


# Both former call sites named this helper `_open_vcf`; keep the name so a
# caller can switch to this module with a single import-line change.
_open_vcf = open_vcf


def drain(buffered) -> None:
    """Consume any input bytes the text/gzip layer left unread.

    Text decoding and gzip framing can stop before the buffered reader has
    hit physical EOF (notably at a gzip member boundary, when a second
    member follows the first). Without this, the hashing stream's digest
    would miss those trailing bytes even though they are part of the file.
    """
    while buffered.read(1 << 20):
        pass


# Canonical (snake_case) genotype-decode failure tags. These are exactly
# `vcf_eligibility`'s current report keys.
INVALID_GENOTYPE = "invalid_genotype"
NON_BIALLELIC_GENOTYPE = "non_biallelic_genotype"
HETEROZYGOUS_GENOTYPE = "heterozygous_genotype"

# `phi_sfs` tags the same three failures with short, hyphenated strings that
# key its own `_GENOTYPE_ERRORS` message table. This maps canonical -> short.
PHI_SFS_REASON_TAGS = {
    INVALID_GENOTYPE: "invalid",
    NON_BIALLELIC_GENOTYPE: "non-biallelic",
    HETEROZYGOUS_GENOTYPE: "heterozygous",
}


def decode_inbred_genotype(gt: str, heterozygous: str) -> int | None | str:
    """Return one inbred individual's allele, decoded from its GT string.

    Each inbred individual contributes a single observed allele, so haploid
    and homozygous diploid calls are accepted and any missing allele makes
    the whole individual uncallable. Returns 0 or 1 for a callable call,
    `None` when the individual is not callable (a missing allele, or a
    heterozygous call under the "missing" policy), or one of
    `INVALID_GENOTYPE`, `NON_BIALLELIC_GENOTYPE`, `HETEROZYGOUS_GENOTYPE`
    when the call itself is malformed rather than merely absent. Callers
    should cache results by GT string, since genotype text is drawn from a
    very small alphabet.
    """
    alleles = gt.replace("|", "/").split("/")
    if not alleles or any(allele == "." for allele in alleles):
        return None
    values: list[int] = []
    for allele in alleles:
        if not allele.isdigit():
            return INVALID_GENOTYPE
        value = int(allele)
        if value not in (0, 1):
            return NON_BIALLELIC_GENOTYPE
        values.append(value)
    if len(set(values)) > 1:
        return None if heterozygous == "missing" else HETEROZYGOUS_GENOTYPE
    return values[0]


def decode_inbred_genotype_short(gt: str, heterozygous: str) -> int | None | str:
    """Same callability decisions as `decode_inbred_genotype`, tagged with
    `phi_sfs`'s short reason strings ("invalid", "non-biallelic",
    "heterozygous") instead of the canonical snake_case ones, for drop-in
    compatibility with `phi_sfs._GENOTYPE_ERRORS`.
    """
    result = decode_inbred_genotype(gt, heterozygous)
    if isinstance(result, str):
        return PHI_SFS_REASON_TAGS[result]
    return result


# Alias under `phi_sfs`'s current private name, so it can import this
# function as `_decode_genotype` with no call-site changes.
_decode_genotype = decode_inbred_genotype_short
