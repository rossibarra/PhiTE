#!/usr/bin/env python3
"""Write the C->T converted reference used for bisulfite mappability.

For every input sequence NAME this writes two records: NAME (forward strand,
C->T) and NAME__rc (reverse complement, C->T).  Reads from either strand of
the original DNA, after bisulfite conversion and read-space C->T collapsing,
are k-mers of this set read in the forward orientation only, so it is indexed
and searched with GenMap's --no-reverse-complement.  N and other non-ACGT
characters are kept as N.

    python mappability/bisulfite_convert.py REF.fa.gz OUT.fa [--chroms chr10 ...]
"""
import argparse
import gzip

FWD = bytes.maketrans(b"ACGTacgtCc", b"ATGTATGTTT")  # upper-case, C->T
COMP = bytes.maketrans(b"ACGTacgt", b"TGCATGCA")
CLEAN = bytes.maketrans(bytes(range(256)), bytes(
    c if chr(c) in "ACGT" else ord("N") for c in range(256)))


def records(path):
    op = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)


def write(out, name, seq):
    out.write(f">{name}\n".encode())
    for i in range(0, len(seq), 100):
        out.write(seq[i:i + 100] + b"\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("out")
    ap.add_argument("--chroms", nargs="*", help="restrict to these sequence names")
    args = ap.parse_args()
    n = 0
    with open(args.out, "wb") as out:
        for name, seq in records(args.fasta):
            if args.chroms and name not in args.chroms:
                continue
            seq = seq.upper().translate(CLEAN)
            write(out, name, seq.translate(FWD))
            write(out, name + "__rc", seq.translate(COMP)[::-1].translate(FWD))
            n += 1
    print(f"converted {n} sequences -> {args.out}")


if __name__ == "__main__":
    main()
