#!/usr/bin/env python3
import re
import sys

if len(sys.argv) != 3:
    print("Usage: python clean_bib.py input.bib output.bib")
    sys.exit(1)

infile  = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r", encoding="utf8", errors="ignore") as f:
    text = f.read()

# -------------------------------------------------------------------
# 1. Normalize various forms of "unknown"
# -------------------------------------------------------------------
patterns_unknown = [
    r'=\s*unknown',
    r'=\s*"unknown"',
    r"=\s*'unknown'",
    r"=\s*''unknown''",
]

for p in patterns_unknown:
    text = re.sub(p, '= {unknown}', text)

# -------------------------------------------------------------------
# 2. Fix empty assignments like: field = ,
# -------------------------------------------------------------------
text = re.sub(r'=\s*,', '= {unknown},', text)

# -------------------------------------------------------------------
# 3. Wrap bare identifiers safely
#    (avoid numeric fields and BibTeX string macros)
# -------------------------------------------------------------------
NUMERIC_FIELDS = {"year", "volume", "number", "pages"}

def wrap_identifier(match):
    field = match.group(1)
    value = match.group(2)

    field_l = field.lower()

    # Skip numeric-only fields (valid bare values)
    if field_l in NUMERIC_FIELDS:
        return match.group(0)

    # Leave already wrapped values alone
    if value.startswith("{") or value.startswith('"'):
        return match.group(0)

    # Do not wrap BibTeX string macros (jan, feb, jfm, etc.)
    if value.isalpha() and value.islower():
        return match.group(0)

    # Wrap simple bare identifiers
    if re.match(r'^[A-Za-z0-9._:-]+$', value):
        return f"{field} = {{{value}}}"

    return match.group(0)

text = re.sub(
    r'(\w+)\s*=\s*([A-Za-z0-9._:-]+)',
    wrap_identifier,
    text
)

# -------------------------------------------------------------------
# 4. Write cleaned output
# -------------------------------------------------------------------
with open(outfile, "w", encoding="utf8") as f:
    f.write(text)

print(f"Cleaned file written to {outfile}")

