#!/usr/bin/env python3
"""
btac_add_doi.py

Normalize a BibTeX file so it is safe for bibtex-autocomplete (btac),
then run btac to add DOI fields using Crossref.

Chunking is OPTIONAL and disabled by default.

Usage:
    python btac_add_doi.py
    python btac_add_doi.py input.bib
    python btac_add_doi.py input.bib output_cleaned.bib
    python btac_add_doi.py input.bib output_cleaned.bib --chunk-size 75

Defaults:
    input.bib          = FDS_general.bib
    output_cleaned.bib = FDS_general_cleaned.bib
    chunking           = OFF

Final output:
    <output_cleaned>.btac.bib

Notes:
- DOI lookup is performed using Crossref only.
- Existing fields are not overwritten.
- '--mark' is used so the script can be safely re-run.
- Chunking may be enabled if rate limiting is encountered on large files.
"""

import argparse
import re
import subprocess
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------
# Cleaning logic
# ---------------------------------------------------------------------
def clean_bib(infile: Path, outfile: Path) -> None:
    with open(infile, "r", encoding="utf8", errors="ignore") as f:
        text = f.read()

    # Normalize various forms of "unknown"
    patterns_unknown = [
        r'=\s*unknown',
        r'=\s*"unknown"',
        r"=\s*'unknown'",
        r"=\s*''unknown''",
    ]
    for p in patterns_unknown:
        text = re.sub(p, '= {unknown}', text)

    # Fix empty assignments like: field = ,
    text = re.sub(r'=\s*,', '= {unknown},', text)

    # Wrap bare identifiers safely
    NUMERIC_FIELDS = {"year", "volume", "number", "pages"}

    def wrap_identifier(match):
        field = match.group(1)
        value = match.group(2)
        field_l = field.lower()

        if field_l in NUMERIC_FIELDS:
            return match.group(0)

        if value.startswith("{") or value.startswith('"'):
            return match.group(0)

        if value.isalpha() and value.islower():
            return match.group(0)

        if re.match(r'^[A-Za-z0-9._:-]+$', value):
            return f"{field} = {{{value}}}"

        return match.group(0)

    text = re.sub(
        r'(\w+)\s*=\s*([A-Za-z0-9._:-]+)',
        wrap_identifier,
        text
    )

    with open(outfile, "w", encoding="utf8") as f:
        f.write(text)


# ---------------------------------------------------------------------
# BibTeX chunk helpers
# ---------------------------------------------------------------------
def split_bib_into_chunks(bibfile: Path, chunk_size: int):
    with open(bibfile, "r", encoding="utf8") as f:
        text = f.read()

    entries = re.split(r'(?=@\w+{)', text)
    entries = [e for e in entries if e.strip()]

    chunks = []
    for i in range(0, len(entries), chunk_size):
        chunk_path = bibfile.with_name(
            f"{bibfile.stem}_chunk_{i//chunk_size:03d}.bib"
        )
        with open(chunk_path, "w", encoding="utf8") as f:
            f.write("".join(entries[i:i + chunk_size]))
        chunks.append(chunk_path)

    return chunks


def merge_bib_files(bib_files, output_file: Path):
    with open(output_file, "w", encoding="utf8") as out:
        for bib in bib_files:
            with open(bib, "r", encoding="utf8") as f:
                out.write(f.read())


# ---------------------------------------------------------------------
# btac invocation
# ---------------------------------------------------------------------
def run_btac(bibfile: Path) -> None:
    """
    Manual equivalent:
        btac --only-complete doi \
             --filter-fields-by-entrytype all \
             --only-query crossref \
             --mark \
             <bibfile>
    """
    cmd = [
        "btac",
        "--only-complete", "doi",
        "--filter-fields-by-entrytype", "all",
        "--only-query", "crossref",
        "--mark",
        str(bibfile),
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Clean a BibTeX file and add DOI fields using btac (Crossref)."
    )
    parser.add_argument(
        "input_bib",
        nargs="?",
        default="FDS_general.bib",
        help="Input BibTeX file (default: FDS_general.bib)",
    )
    parser.add_argument(
        "output_cleaned",
        nargs="?",
        default="FDS_general_cleaned.bib",
        help="Cleaned BibTeX file (default: FDS_general_cleaned.bib)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=0,
        help="Optional chunk size for btac processing (default: no chunking)",
    )
    parser.add_argument(
        "--sleep",
        type=int,
        default=45,
        help="Sleep time between chunks in seconds (default: 45)",
    )

    args = parser.parse_args()

    input_bib = Path(args.input_bib)
    output_cleaned = Path(args.output_cleaned)

    if not input_bib.exists():
        print(f"ERROR: input file not found: {input_bib}", file=sys.stderr)
        sys.exit(1)

    print(f"Cleaning {input_bib} → {output_cleaned}")
    clean_bib(input_bib, output_cleaned)

    final_btac = output_cleaned.with_suffix(".btac.bib")

    if args.chunk_size > 0:
        print(f"\nRunning btac with chunking (size={args.chunk_size})")
        chunks = split_bib_into_chunks(output_cleaned, args.chunk_size)
        outputs = []

        for i, chunk in enumerate(chunks, start=1):
            print(f"\n[{i}/{len(chunks)}] btac on {chunk.name}")
            run_btac(chunk)
            outputs.append(chunk.with_suffix(".btac.bib"))

            if i < len(chunks):
                print(f"Sleeping {args.sleep} s...")
                time.sleep(args.sleep)

        merge_bib_files(outputs, final_btac)

    else:
        print("\nRunning btac without chunking")
        run_btac(output_cleaned)

    print("\nDone.")
    print(f"Cleaned file:      {output_cleaned}")
    print(f"DOI-enriched file: {final_btac}")


if __name__ == "__main__":
    main()
