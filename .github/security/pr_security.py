#!/usr/bin/env python3
"""Trusted PR admission and ClamAV helpers. Never import or execute PR code."""

import argparse
import base64
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile
import tempfile
import time
import urllib.error
import urllib.request


ORGANIZATION = "firemodels"
ADMISSION = "FDS / PR admission"
SCAN = "FDS / ClamAV scan"
BUILDS = "FDS / PR builds"
MAX_DOWNLOAD = 1024 * 1024 * 1024
MAX_FILE = 512 * 1024 * 1024
MAX_TOTAL = 4 * 1024 * 1024 * 1024
MAX_FILES = 100000
LFS_POINTER = b"version https://git-lfs.github.com/spec/v1"
SHA = re.compile(r"[0-9a-f]{40}\Z")
REPOSITORY = re.compile(r"[A-Za-z0-9][A-Za-z0-9-]{0,38}/[A-Za-z0-9_.-]+\Z")


def api(path, data=None, authenticated=True):
    if not path.startswith("/"):
        raise ValueError("Expected a GitHub API path")
    headers = {"Accept": "application/vnd.github+json", "User-Agent": "fds-pr-security",
               "X-GitHub-Api-Version": "2022-11-28"}
    if authenticated:
        headers["Authorization"] = "Bearer " + os.environ["GH_TOKEN"]
    body = None if data is None else json.dumps(data).encode()
    request = urllib.request.Request("https://api.github.com" + path, data=body, headers=headers)
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read()
        return json.loads(content) if content else None


def organization_member(user):
    # The public membership endpoint works in personal forks without an org secret.
    # Private membership, unavailable APIs, and bots do not grant a fast track.
    login = user.get("login", "")
    if user.get("type") != "User" or not re.fullmatch(r"[A-Za-z0-9-]{1,39}", login):
        return False
    try:
        return api(f"/orgs/{ORGANIZATION}/public_members/{login}", authenticated=False) is None
    except (urllib.error.URLError, TimeoutError, ValueError):
        print(f"Public membership could not be verified for {login}; scanning is required.")
        return False


def validate_revision(repository, sha):
    if (not REPOSITORY.fullmatch(repository) or repository.rsplit("/", 1)[-1] in (".", "..")
            or not SHA.fullmatch(sha)):
        raise ValueError("Invalid repository or commit SHA")


def status(repository, sha, context, state, description):
    validate_revision(repository, sha)
    api(f"/repos/{repository}/statuses/{sha}", {
        "context": context, "state": state, "description": description,
        "target_url": f"https://github.com/{repository}/actions/runs/{os.environ['GITHUB_RUN_ID']}",
    })


def outputs(values):
    with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as stream:
        for name, value in values.items():
            value = str(value)
            if "\n" in value or "\r" in value:
                raise ValueError("Invalid workflow output")
            stream.write(f"{name}={value}\n")


def build_selection(files, complete):
    paths = [name for item in files for name in (item["filename"], item.get("previous_filename", ""))]
    conventional = not complete or any(
        name.startswith((".github/", "Build/", "Source/", "Utilities/Python/")) for name in paths
    )
    cmake = conventional or any(name in ("CMakeLists.txt", "CMakePresets.json") for name in paths)
    return conventional, cmake


def prepare():
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
    repository = os.environ["GITHUB_REPOSITORY"]
    number = int(event["number"])
    head_sha = event["pull_request"]["head"]["sha"]
    validate_revision(repository, head_sha)
    status(repository, head_sha, ADMISSION, "pending", "Checking identity and the exact PR revision")
    status(repository, head_sha, SCAN, "pending", "ClamAV scan is pending")
    status(repository, head_sha, BUILDS, "pending", "Builds are waiting for admission")
    try:
        pull = api(f"/repos/{repository}/pulls/{number}")
        if pull["state"] != "open" or pull["head"]["sha"] != head_sha:
            raise ValueError("The PR changed or closed; this event is stale")
        head_repo = pull["head"]["repo"]["full_name"]
        validate_revision(head_repo, head_sha)
        if pull["head"]["repo"]["private"]:
            raise ValueError("This scanner requires publicly downloadable source snapshots")

        # Never use author_association, fork ownership, commit author strings, or
        # triggering_actor (the person re-running a job) as a trust decision.
        trusted = organization_member(pull["user"]) and organization_member(event["sender"])
        base_sha = pull["base"]["sha"]
        merge_sha = ""
        for attempt in range(10):
            candidate = pull.get("merge_commit_sha") or ""
            if SHA.fullmatch(candidate):
                commit = api(f"/repos/{repository}/git/commits/{candidate}")
                parents = [parent["sha"] for parent in commit["parents"]]
                if parents == [base_sha, head_sha]:
                    merge_sha = candidate
                    break
            if pull.get("mergeable") is False:
                break
            time.sleep(3)
            pull = api(f"/repos/{repository}/pulls/{number}")
            if pull["state"] != "open" or pull["head"]["sha"] != head_sha or pull["base"]["sha"] != base_sha:
                raise ValueError("The PR changed while preparing the scan")

        files = []
        for page in range(1, 31):
            batch = api(f"/repos/{repository}/pulls/{number}/files?per_page=100&page={page}")
            files.extend(batch)
            if len(batch) < 100:
                break
        conventional, cmake = build_selection(files, len(files) == pull["changed_files"])
        outputs({"head_sha": head_sha, "head_repo": head_repo, "base_sha": base_sha,
                 "merge_sha": merge_sha, "trusted": str(trusted).lower(),
                 "conventional": str(conventional).lower(), "cmake": str(cmake).lower()})
        print(f"PR #{number}: head={head_sha}, merge={merge_sha or 'unavailable'}, trusted={trusted}")
    except Exception:
        status(repository, head_sha, ADMISSION, "failure", "Admission preparation failed; review workflow logs")
        status(repository, head_sha, SCAN, "error", "Scan could not be prepared")
        status(repository, head_sha, BUILDS, "error", "Builds could not be prepared")
        raise


def download(repository, sha, destination):
    validate_revision(repository, sha)
    # Public codeload requests carry no GitHub token and execute no checkout hooks.
    url = f"https://codeload.github.com/{repository}/tar.gz/{sha}"
    total = 0
    digest = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=120) as response, destination.open("wb") as stream:
        while chunk := response.read(1024 * 1024):
            total += len(chunk)
            if total > MAX_DOWNLOAD:
                raise ValueError("Source archive exceeds the download limit; scan is incomplete")
            digest.update(chunk)
            stream.write(chunk)
    print(f"Snapshot {repository}@{sha}, archive SHA256={digest.hexdigest()}")


def snapshot_manifest(repository, sha):
    """Get the full Git tree, independent of archive export/LFS settings."""
    validate_revision(repository, sha)
    commit = api(f"/repos/{repository}/git/commits/{sha}")
    tree_sha = commit["tree"]["sha"]
    validate_revision(repository, tree_sha)
    tree = api(f"/repos/{repository}/git/trees/{tree_sha}?recursive=1")
    if tree.get("truncated", True) or tree["sha"] != tree_sha:
        raise ValueError("Incomplete Git tree; cannot verify scan coverage")
    files = {}
    total = 0
    for entry in tree["tree"]:
        if entry["type"] == "tree" and entry["mode"] == "040000":
            continue
        if entry["type"] != "blob" or entry["mode"] not in ("100644", "100755"):
            raise ValueError("Git links, submodules, or special files require manual review")
        name = entry["path"]
        size = entry["size"]
        if name in files or size < 0 or size > MAX_FILE or not SHA.fullmatch(entry["sha"]):
            raise ValueError("Invalid or oversized Git tree entry")
        files[name] = (size, entry["sha"])
        total += size
        if len(files) > MAX_FILES or total > MAX_TOTAL:
            raise ValueError("Git tree exceeds scanning limits")
    if not files:
        raise ValueError("Empty Git tree")
    return files


def verify_snapshot(destination, manifest):
    files = {p.relative_to(destination).as_posix(): p for p in destination.rglob("*") if p.is_file()}
    if files.keys() != manifest.keys():
        raise ValueError("Archive omits or adds Git files (for example export-ignore); scan is incomplete")
    for name, path in files.items():
        size, expected_hash = manifest[name]
        if path.stat().st_size != size:
            raise ValueError("Archive content differs from the Git blob")
        digest = hashlib.sha1(f"blob {size}\0".encode())
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
        if digest.hexdigest() != expected_hash:
            raise ValueError("Archive content differs from the Git blob (for example export-subst or LFS)")
    print(f"Verified all {len(files)} files against Git blob hashes")


def extract_snapshot(archive, destination):
    """Extract only regular files, without restoring permissions or following links."""
    total = 0
    count = 0
    root = None
    seen = set()
    with tarfile.open(archive, "r|gz") as source:
        for member in source:
            path = PurePosixPath(member.name)
            if (path.is_absolute() or ".." in path.parts or not path.parts
                    or "\\" in member.name or any(ord(char) < 32 or ord(char) == 127 for char in member.name)):
                raise ValueError("Unsafe archive path; scan is incomplete")
            root = root or path.parts[0]
            if path.parts[0] != root:
                raise ValueError("Unexpected archive root")
            if member.isdir():
                continue
            if not member.isfile() or len(path.parts) < 2:
                raise ValueError("Links or special files need manual review; scan is incomplete")
            relative = Path(*path.parts[1:])
            if relative in seen:
                raise ValueError("Duplicate archive path")
            seen.add(relative)
            count += 1
            total += member.size
            if member.size < 0 or member.size > MAX_FILE or total > MAX_TOTAL or count > MAX_FILES:
                raise ValueError("Source exceeds scanning limits; scan is incomplete")
            target = destination / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            with source.extractfile(member) as incoming, target.open("xb") as outgoing:
                prefix = incoming.read(min(member.size, 1024))
                if prefix.startswith(LFS_POINTER):
                    raise ValueError("Unresolved Git LFS object; scan is incomplete")
                if relative.name == ".gitmodules" and prefix.strip():
                    raise ValueError("Submodule content needs a separate scan; scan is incomplete")
                outgoing.write(prefix)
                written = len(prefix)
                while chunk := incoming.read(1024 * 1024):
                    written += len(chunk)
                    if written > member.size:
                        raise ValueError("Archive file size mismatch")
                    outgoing.write(chunk)
                if written != member.size:
                    raise ValueError("Truncated archive file")
    if count == 0:
        raise ValueError("Empty source snapshot")
    print(f"Prepared {count} files ({total} bytes) for scanning")
    return count


def clamav(path):
    expected = sum(item.stat().st_size > 0 for item in path.rglob("*") if item.is_file())
    result = subprocess.run([
        "clamscan", "--recursive=yes", "--infected", "--alert-exceeds-max=yes",
        "--alert-encrypted=yes", "--max-filesize=512M", "--max-scansize=2048M",
        "--max-recursion=40", "--max-files=100000", "--max-scantime=120000",
        "--follow-dir-symlinks=0", "--follow-file-symlinks=0", str(path),
    ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1200)
    # Prefix untrusted filenames so they cannot be interpreted as workflow commands.
    for line in result.stdout.splitlines():
        print("ClamAV | " + line)
    if result.returncode != 0 or re.search(r"\b(?:ERROR|WARNING)\b", result.stdout):
        raise RuntimeError(f"ClamAV did not complete cleanly (exit {result.returncode})")
    match = re.search(r"^Scanned files:\s+(\d+)\s*$", result.stdout, re.MULTILINE)
    if match is None or int(match[1]) == 0 or int(match[1]) < expected:
        raise RuntimeError("ClamAV did not report covering all nonempty source files")


def scanner_self_test(directory):
    # Generate the standard harmless EICAR test file only on the disposable runner.
    # Keep the test bytes out of normal source checkouts and workstation test runs.
    sample = base64.b64decode(
        "WDVPIVAlQEFQWzRcUFpYNTQoUF4pN0NDKTd9JEVJQ0FSLVNUQU5EQVJELUFOVElWSVJVUy1URVNULUZJTEUhJEgrSCo="
    )
    probe = directory / "scanner-probe"
    probe.write_bytes(sample)
    result = subprocess.run(["clamscan", "--infected", str(probe)], capture_output=True, text=True, timeout=120)
    probe.unlink()
    if result.returncode != 1 or "eicar" not in result.stdout.lower():
        raise RuntimeError("ClamAV failed the EICAR detection self-test")
    print("ClamAV detection self-test passed")


def scan():
    with tempfile.TemporaryDirectory(prefix="fds-virus-scan-") as name:
        work = Path(name)
        scanner_self_test(work)
        snapshots = [(os.environ["HEAD_REPO"], os.environ["HEAD_SHA"])]
        merge_sha = os.environ.get("MERGE_SHA", "")
        if merge_sha:
            snapshots.append((os.environ["GITHUB_REPOSITORY"], merge_sha))
        for index, (repository, sha) in enumerate(snapshots):
            archive = work / f"snapshot-{index}.tar.gz"
            destination = work / f"snapshot-{index}"
            destination.mkdir()
            manifest = snapshot_manifest(repository, sha)
            download(repository, sha, archive)
            extract_snapshot(archive, destination)
            verify_snapshot(destination, manifest)
            clamav(destination)
        if not merge_sha:
            raise RuntimeError("PR head scanned, but no current merge snapshot is available; admission remains blocked")


def report(kind):
    repository = os.environ["GITHUB_REPOSITORY"]
    sha = os.environ["HEAD_SHA"]
    trusted = os.environ.get("TRUSTED") == "true"
    external = os.environ.get("EXTERNAL_RESULT", "")
    if kind == "admission":
        accepted = bool(os.environ["MERGE_SHA"]) and (trusted or external == "success")
        if accepted:
            number = int(os.environ["PR_NUMBER"])
            current = api(f"/repos/{repository}/pulls/{number}")
            if (current["state"] != "open" or current["head"]["sha"] != sha
                    or current["base"]["sha"] != os.environ["BASE_SHA"]):
                # Do not let an older run replace a newer run's pending status.
                raise RuntimeError("PR changed before admission; a new scan is required")
        context = ADMISSION
        description = ("Verified firemodels member; scan continues in background" if trusted else
                       "Complete PR head and merge snapshot scans passed") if accepted else "PR not admitted; scan or merge snapshot is incomplete"
    elif kind == "builds":
        results = json.loads(os.environ["BUILD_RESULTS"])
        conventional = os.environ["CONVENTIONAL"] == "true"
        cmake = os.environ["CMAKE"] == "true"
        expected = {"prepare": True, "admission": True, "line-endings": True,
                    "linux": conventional, "macos": conventional, "windows": conventional, "cmake": cmake}
        accepted = all(results.get(name, {}).get("result") == ("success" if selected else "skipped")
                       for name, selected in expected.items())
        context = BUILDS
        description = "All required PR builds/checks passed" if accepted else "PR builds/checks failed or were blocked"
    else:
        accepted = (os.environ["TRUSTED_RESULT"] if trusted else external) == "success"
        context = SCAN
        description = "Complete PR head and merge snapshot scans passed" if accepted else "ClamAV scan failed or was incomplete; inspect workflow logs"
    status(repository, sha, context, "success" if accepted else "failure", description)
    if not accepted:
        raise RuntimeError(description)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "scan", "admission", "report", "builds"))
    command = parser.parse_args().command
    if command == "prepare":
        prepare()
    elif command == "scan":
        scan()
    else:
        report(command)


if __name__ == "__main__":
    main()
