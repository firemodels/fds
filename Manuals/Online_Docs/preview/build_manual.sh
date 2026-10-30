#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Honor an explicit interpreter; otherwise prefer the project environment.
if [[ -n "${FDS_DOCS_PYTHON:-}" ]]; then
  exec "$FDS_DOCS_PYTHON" "$script_dir/build_preview.py" --full "$@"
fi
if [[ -n "${VIRTUAL_ENV:-}" && -x "$VIRTUAL_ENV/bin/python3" ]]; then
  if ! "$VIRTUAL_ENV/bin/python3" -c 'import lxml' >/dev/null 2>&1; then
    printf '%s\n' 'The active Python environment needs lxml. Run: python -m pip install lxml==6.1.1' >&2
    exit 1
  fi
  exec "$VIRTUAL_ENV/bin/python3" "$script_dir/build_preview.py" --full "$@"
fi
for docs_python in "$script_dir/../.venv/bin/python3" python3 \
  "$HOME/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3"; do
  if "$docs_python" -c 'import lxml' >/dev/null 2>&1; then
    exec "$docs_python" "$script_dir/build_preview.py" --full "$@"
  fi
done
printf '%s\n' 'Python with lxml is required. Install requirements.txt from the repository .github directory in a virtual environment, or set FDS_DOCS_PYTHON.' >&2
exit 1
