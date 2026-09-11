# PR virus scanning and admission

`pr-security.yml` runs for every pull request against this repository, including
pull requests against a developer's fork. The receiving repository must have the
workflow and helpers on its trusted base branch and have GitHub Actions enabled.

The gate checks the PR author and the event sender against the **public membership
of the `firemodels` GitHub organization**. Both must be verified members for the
fast path. Fork ownership, repository collaborator status, previous contributions,
commit author text, and the person re-running a workflow do not grant this trust.
Private memberships cannot be verified without additional organization access;
private members and API failures therefore take the scanning-required path.
No organization secret or personal access token is needed.

| Contributor | `FDS / PR admission` | `FDS / ClamAV scan` |
| --- | --- | --- |
| Verified organization member | Passes after preparing a current merge snapshot; builds may start immediately | Runs independently in the background |
| Everyone else | Passes only after both complete source snapshots scan successfully | Must finish successfully before builds can start |

The snapshots are the exact submitted head commit and the proposed merge commit.
The gate verifies the merge commit's parent SHAs. Builds check out that same merge
SHA, rather than a moving branch or PR ref. Updates start a new run and cancel the
old one. A missing merge snapshot, conflict, stale event, or scanner error does
not admit outside code. A conflicted PR's head is still scanned, but the combined
scan cannot pass until a current merge snapshot is available.

The existing Linux, macOS, Windows, CMake, and line-ending workflows are callable
by the gate instead of running directly on `pull_request`. Their push behavior is
preserved. PR source is never executed in the identity or scanning jobs. Only the
preparation and status-reporting jobs receive `statuses: write`; scan and build
jobs have `contents: read`, do not inherit secrets, and do not retain checkout
credentials. The workflow/helper checkouts in the security jobs are pinned to
the trusted workflow revision, not to the PR's files.

## Required setup in upstream AND every developer fork

These repository settings are essential and are **not installed by copying YAML**:

1. Enable GitHub Actions in the receiving repository and install these files on
   each branch that accepts PRs. Until that is done, its existing workflows and
   merge rules still apply. The workflow deliberately follows the local trusted
   revision so forks do not depend on a moving upstream workflow.
2. Protect the target branches using a branch rule or ruleset. Require the exact
   commit status **`FDS / PR admission`**, require branches to be up to date before
   merging, and restrict bypass/direct pushes as appropriate. Also retain the
   project's review requirements. Use **`FDS / PR builds`** for the existing
   platform builds and line-ending checks: it reports their aggregate result on
   the PR head. Replace old required job names from `pull_request` workflows;
   `pull_request_target` job checks attach to the base revision instead. The build
   status also does not wait for member scans. Do **not** require
   `FDS / ClamAV scan` if members must be able to merge before their scan finishes.
3. In Actions settings, require approval for fork-PR workflows from **all outside
   collaborators**, not just first-time contributors. Do not approve PR-supplied
   workflows until admission succeeds and any workflow changes have been reviewed.
   An attacker can add a *new* `pull_request` workflow in their PR; the trusted gate
   cannot stop GitHub from scheduling unrelated workflows through YAML alone.
4. Protect edits to `.github/workflows/` and `.github/security/` with maintainer
   review. Where your GitHub plan supports organization-required workflows, use a
   ruleset requiring the trusted admission workflow as an additional control.
   A status name alone is not cryptographic proof of which workflow produced it;
   do not allow arbitrary PR workflows to publish trusted statuses.
5. Keep unreviewed PRs off self-hosted/NIST runners. This workflow uses disposable
   standard GitHub-hosted runners. Separate external bots or developer machines
   must independently check admission for the exact head SHA before downloading
   and executing a PR. Updating these workflows does not gate an external Firebot.

For merge-queue deployments, add and test an equivalent scan for the queue's
`merge_group` snapshot before enabling the queue. This configuration covers ordinary
PR merges, not merge queues.

## Scan scope and limits

ClamAV is installed from the Ubuntu package repositories, and `freshclam` must
successfully update its signatures on every run. A harmless EICAR self-test must
be detected before candidate files are scanned. The self-test is generated only
inside the disposable runner; it is not committed as a test file.

The scanner downloads public GitHub source archives without credentials and scans
all regular files, including hidden files, documentation, scripts, and committed
binaries. It never sources a script or builds the submitted code. Archives are
extracted without restoring executable bits or following links. Both the scan
summary and the downloaded archive SHA256 are retained in the Actions logs.
The read-only token retrieves the complete Git tree, and every extracted file
must match its Git blob hash. Truncated trees, files omitted by `export-ignore`,
content changed by `export-subst`, and gitlinks without a `.gitmodules` file
cannot silently pass as a complete scan. Source archive downloads carry no token.

The scanner rejects path traversal, symbolic/hard links, special files, unresolved
Git LFS pointers, and submodules instead of reporting an incomplete scan as clean.
Limits are 1 GiB per downloaded archive, 512 MiB per source file, 4 GiB of extracted
data, and 100,000 source files. ClamAV's nested-content limits, encrypted-content
alerts, warnings, nonzero exit statuses, missing file coverage, and timeouts fail
the scan. Exceeding a limit requires a reviewed change to the scanner or a separate
security review; do not bypass the gate by excluding the file.

ClamAV detects recognizable malware; a successful scan does not prove that code
is benign or that scientific results are correct. Builds can still execute harmful
code that antivirus misses. Code review and isolated build environments remain
necessary. A compromised organization-member account also retains its fast path.

## Validation and rollout

Run the offline tests from the repository root:

```sh
python3 -m unittest discover -s .github/security -p 'test_*.py' -v
```

These tests use mocked GitHub/ClamAV calls, synthetic non-malicious archives, and
Ruby's YAML parser. They do not connect to GitHub, scan with a real engine, or
write an EICAR file on a workstation. They also run in the existing Linux Python
job after admission.

Before relying on the gate, validate it in a developer fork on GitHub:

- A clean organization-member PR admits promptly while its scan continues.
- A clean nonmember PR starts no build jobs until its scan succeeds.
- The runner's EICAR self-test passes, and a controlled EICAR PR from a nonmember
  blocks admission and all build jobs. Use an isolated test repository for this.
- A scanner/signature-update failure blocks nonmembers; member scans still report
  their failure separately.
- A new PR commit requires a fresh status, and repository rules prevent merging
  a pending or failed nonmember PR.

Do not use production secrets or self-hosted runners during these rollout tests.
