# PR author policy and virus scanning

The policy is based on the **PR author's** numeric GitHub account ID in
`trusted-developers.json`. The person pushing, editing, approving, or rerunning
the PR does not change its classification. Usernames in the list are labels;
renaming an account does not change its numeric ID. Missing or null IDs are not
trusted. The list is always read from the receiving repository's trusted revision.

| PR author | Merge permission | Builds and tests | Virus scan |
| --- | --- | --- | --- |
| Listed developer | Admission passes immediately, subject to other branch rules | May start immediately | Continues independently, including after merge |
| Anyone else | Admission stays pending until the scan passes | Wait for successful admission | Failure, cancellation, or incomplete scan blocks admission |

## Two workflows with separate responsibilities

**`pr-security.yml`** runs the trusted policy under `pull_request_target`:

1. **Check PR author:** publish `FDS / PR admission` immediately for a listed
   author, or leave it pending for everyone else.
2. **Virus scan:** scan the exact PR head and proposed merge as data. This job has
   a read-only token and never executes PR scripts.
3. **Report virus scan:** publish `FDS / ClamAV scan`. For an unlisted author, also
   pass or fail admission. A developer's admission is not changed by the scan.

Only the author and reporting jobs can write commit statuses. All three jobs
check out the trusted workflow revision, never PR code. Their narrowly scoped
`allow-unsafe-checkout` setting accommodates checkout v7's restriction on fork
PRs in `pull_request_target`. This setting is absent from every build job.

**`pr-builds.yml`** runs under ordinary `pull_request`, with read-only permissions:

1. Wait for admission of the **exact merge commit** supplied by GitHub. The wait
   runs only the trusted admission reader from the base commit, before checking
   out any PR code. It checks that admission references the security workflow.
2. Run the selected Linux, macOS, Windows, CMake, and line-ending checks. The
   existing push workflows remain available. PR builds use normal PR cache scope.
3. Report the ordinary check `FDS / PR builds`; no write token is needed.

Admission is published on the head commit for merge rules and on the merge
commit for build coordination. A changed head or base cannot reuse a result for
a different merge. The wait fails after 50 minutes if admission never arrives;
inspect the security run and rerun builds if necessary. Listed authors wait only
for the short author check, not the scanner. A scan still reports its result if
the PR was merged while it ran.

## Repository settings and rollout

YAML alone cannot enforce merge rules or stop an unreviewed PR from rewriting
its own `pull_request` build workflow. Configure these settings in upstream and
in each developer fork that accepts PRs:

- Require **`FDS / PR admission`** in the target branch's protection/ruleset, and
  require the branch to be up to date. Restrict bypass/direct pushes as needed.
- Keep existing build/review requirements. If builds are required, select
  **`FDS / PR builds`**, now an ordinary check from the PR builds workflow instead
  of a commit status published by the old gate. Verify the rule after rollout.
- **Do not require `FDS / ClamAV scan` or the whole PR security workflow** if
  listed developers must be able to merge while scanning is still running.
- Retain GitHub's workflow approval requirement for **all outside contributors**.
  Before approving a workflow changed by an outside PR, wait for admission and
  review its workflow changes. The scan workflow itself uses trusted base code.
  This GitHub approval is separate from the developer list and may still apply
  to a listed author who is not a repository collaborator.
- Protect edits to `.github/workflows/` and `.github/security/` with maintainer
  review. Keep PR jobs on disposable GitHub-hosted runners with read-only tokens
  and no inherited secrets. Unrelated external bots such as Firebot need their
  own admission check; this policy only gates these GitHub Actions builds.

Install the revised policy, helper, and build workflow together on the receiving
base branch. A PR's edited policy does not govern its own admission, and rerunning
an old workflow revision does not install the update. After installation, trigger
a new PR event and verify the required checks in branch protection. The old
`virus-scan.yml` reusable workflow and the gate's build-dispatch/reporting jobs
have been removed.

This covers ordinary PRs, not merge queues. A scan is not a guarantee that code
is safe; approval and review remain necessary. A compromised listed developer
account still receives the fast path. To add an author, obtain their numeric ID
with `gh api users/LOGIN --jq '{login, id}'` and review the list change normally.

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
necessary. A compromised listed developer account also retains its fast path.

## Validation

From the repository root:

```sh
python3 -m unittest discover -s .github/security -p 'test_*.py' -v
```

The tests mock GitHub and ClamAV; they do not create malware, access the network,
or execute submitted PR files. They cover author-only trust, early developer
admission, failed scans, stale revisions, exact-merge build admission, permissions,
and required build results. The Linux Python job also runs them.

After rollout, check a listed-author PR and an outside-author PR on GitHub. For
the first, builds/admission should pass before scanning finishes; for the second,
builds must wait. Confirm a failed scanner blocks outside admission while leaving
developer admission intact. These live checks require GitHub access and are not
replaced by the offline tests.
