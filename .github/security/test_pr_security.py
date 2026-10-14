"""Offline regression tests; no real malware/EICAR files or network requests."""

import copy
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch
import urllib.error

import pr_security as security


HEAD = "a" * 40
BASE = "b" * 40
MERGE = "c" * 40
ROOT = Path(__file__).resolve().parents[2]


def pull_request():
    return {"state": "open", "head": {"sha": HEAD, "repo": {"full_name": "outsider/fds", "private": False}},
            "base": {"sha": BASE}, "user": {"login": "author", "type": "User"},
            "merge_commit_sha": MERGE, "mergeable": True, "changed_files": 1}


class IdentityTests(unittest.TestCase):
    def test_matching_id_retains_trust_after_rename_but_reused_username_does_not(self):
        with patch.object(security, "api") as api:
            for login, account_id, expected in (("member", 123, True), ("renamed", 123, True), ("member", 456, False)):
                user = {"login": login, "id": account_id, "type": "User"}
                self.assertEqual(security.trusted_developer(user, {123}), expected)
            api.assert_not_called()

    def test_bots_and_missing_or_invalid_ids_are_not_trusted(self):
        for account_id in (None, "123", True, 456):
            self.assertFalse(security.trusted_developer({"id": account_id, "type": "User"}, {123}))
        self.assertFalse(security.trusted_developer({"id": 123, "type": "Bot"}, {123}))

    def test_null_ids_are_excluded_and_invalid_configuration_fails(self):
        with tempfile.TemporaryDirectory() as name:
            path = Path(name) / "trusted-developers.json"
            with patch.object(security, "TRUSTED_DEVELOPERS", path):
                path.write_text(json.dumps({"developers": [{"id": 123}, {"id": None}]}))
                self.assertEqual(security.load_developers(), {123})
                for account_id in (True, "123", 0, -1):
                    path.write_text(json.dumps({"developers": [{"id": account_id}]}))
                    with self.assertRaises(ValueError):
                        security.load_developers()

    def test_policy_is_loaded_beside_trusted_helper_not_from_working_directory(self):
        self.assertEqual(security.TRUSTED_DEVELOPERS, Path(security.__file__).with_name("trusted-developers.json"))
        with tempfile.TemporaryDirectory() as name:
            original = Path.cwd()
            expected = security.load_developers()
            try:
                os.chdir(name)
                Path("trusted-developers.json").write_text('{"developers": [{"login": "attacker", "id": 999}]}')
                self.assertEqual(security.load_developers(), expected)
            finally:
                os.chdir(original)


class PreparationTests(unittest.TestCase):
    def prepare(self, pull=None, memberships=(False, False), event_sha=HEAD, parents=None, later_pull=None):
        pull = pull or pull_request()
        sequence = [pull, later_pull] if later_pull is not None else [pull]
        event = {"number": 7, "pull_request": {"head": {"sha": event_sha}},
                 "sender": {"login": "sender", "type": "User"}}

        def api(path, **kwargs):
            if path.endswith("/git/commits/" + MERGE):
                return {"parents": [{"sha": sha} for sha in (parents if parents is not None else [BASE, HEAD])]}
            if "/files?" in path:
                return [{"filename": "Source/main.f90"}]
            return copy.deepcopy(sequence.pop(0) if len(sequence) > 1 else sequence[0])

        with tempfile.TemporaryDirectory() as folder:
            event_path = Path(folder) / "event.json"
            event_path.write_text(json.dumps(event))
            with patch.dict(os.environ, {"GITHUB_EVENT_PATH": str(event_path), "GITHUB_REPOSITORY": "developer/fds"}), \
                    patch.object(security, "api", side_effect=api), \
                    patch.object(security, "status") as status, \
                    patch.object(security, "trusted_developer", side_effect=memberships), \
                    patch.object(security, "outputs") as output, patch.object(security.time, "sleep"):
                security.prepare()
                self.preparation_statuses = status.call_args_list
                self.assertEqual(status.call_args_list[0].args[:2], ("developer/fds", event_sha))
                return output.call_args.args[0]

    def test_fork_owner_or_existing_collaborator_is_not_automatically_trusted(self):
        pull = pull_request()
        pull["author_association"] = "OWNER"
        result = self.prepare(pull)
        self.assertEqual(result["trusted"], "false")
        self.assertEqual(result["merge_sha"], MERGE)

    def test_only_pr_author_determines_trust(self):
        self.assertEqual(self.prepare(memberships=(True, True))["trusted"], "true")
        self.assertEqual(self.prepare(memberships=(True, False))["trusted"], "true")
        self.assertEqual(self.prepare(memberships=(False, True))["trusted"], "false")

    def test_developer_gets_head_and_merge_admission_before_scanning(self):
        self.prepare(memberships=(True, False))
        successes = [(call.args[1], call.args[2]) for call in self.preparation_statuses
                     if call.args[3] == "success"]
        self.assertEqual(successes, [(MERGE, security.ADMISSION), (HEAD, security.ADMISSION)])
        scans = [call for call in self.preparation_statuses if call.args[2] == security.SCAN]
        self.assertEqual([call.args[3] for call in scans], ["pending"])

    def test_outsider_has_no_admission_until_scanner_finishes(self):
        self.prepare(memberships=(False, True))
        self.assertTrue(all(call.args[3] == "pending" for call in self.preparation_statuses))
        admissions = [call.args[1] for call in self.preparation_statuses if call.args[2] == security.ADMISSION]
        self.assertEqual(admissions, [HEAD, MERGE])

    def test_stale_event_does_not_admit_a_new_revision(self):
        with self.assertRaisesRegex(ValueError, "stale"):
            self.prepare(event_sha="d" * 40)

    def test_merged_pr_can_prepare_scan_without_resetting_gate_statuses(self):
        pull = pull_request()
        pull.update(state="closed", merged=True)
        for parents in ([BASE, HEAD], [BASE]):
            with self.subTest(parents=parents):
                result = self.prepare(pull=pull, memberships=(True, True), parents=parents)
                self.assertEqual(result["head_sha"], HEAD)
                self.assertEqual(result["merge_sha"], MERGE)
                self.assertEqual([call.args[2] for call in self.preparation_statuses], [security.SCAN])

    def test_merge_during_preparation_keeps_scanning(self):
        initial = pull_request()
        initial["merge_commit_sha"] = None
        merged = pull_request()
        merged.update(state="closed", merged=True)
        merged["base"]["sha"] = MERGE
        result = self.prepare(pull=initial, later_pull=merged, parents=[BASE])
        self.assertEqual(result["merge_sha"], MERGE)
        self.assertFalse(any(call.args[3] in ("failure", "error") for call in self.preparation_statuses))

    def test_open_pr_still_requires_matching_merge_parents(self):
        self.assertEqual(self.prepare(parents=[BASE])["merge_sha"], "")

    def test_merged_pr_with_deleted_source_fork_uses_receiving_repository(self):
        pull = pull_request()
        pull.update(state="closed", merged=True)
        pull["head"]["repo"] = None
        self.assertEqual(self.prepare(pull=pull)["head_repo"], "developer/fds")

    def test_conflict_still_permits_head_scan_but_no_build_snapshot(self):
        pull = pull_request()
        pull.update(merge_commit_sha=None, mergeable=False)
        self.assertEqual(self.prepare(pull)["merge_sha"], "")

    def test_path_filters_include_python_and_renames_and_fail_open_to_more_tests(self):
        self.assertEqual(security.build_selection([{"filename": "Manuals/a.tex"}], True), (False, False))
        self.assertEqual(security.build_selection([{"filename": "Utilities/Python/setup.py"}], True), (True, True))
        self.assertEqual(security.build_selection([{"filename": "CMakeLists.txt"}], True), (False, True))
        self.assertEqual(security.build_selection([{"filename": "old.txt", "previous_filename": "Build/makefile"}], True), (True, True))
        self.assertEqual(security.build_selection([], False), (True, True))


class AdmissionTests(unittest.TestCase):
    def finish(self, trusted=False, outcome="success", current=None, merge=MERGE):
        environment = {"GITHUB_REPOSITORY": "developer/fds", "HEAD_SHA": HEAD, "BASE_SHA": BASE,
                       "PR_NUMBER": "7", "MERGE_SHA": merge, "TRUSTED": str(trusted).lower(),
                       "SCAN_RESULT": outcome}
        with patch.dict(os.environ, environment), patch.object(security, "api", return_value=current or pull_request()), \
                patch.object(security, "status") as status:
            try:
                security.finish()
                passed = True
            except RuntimeError:
                passed = False
            return passed, status.call_args_list

    def test_outsiders_require_successful_scan_on_head_and_exact_merge(self):
        passed, calls = self.finish()
        self.assertTrue(passed)
        self.assertEqual([(call.args[1], call.args[2], call.args[3]) for call in calls],
                         [(HEAD, security.SCAN, "success"), (MERGE, security.ADMISSION, "success"),
                          (HEAD, security.ADMISSION, "success")])
        for result in ("failure", "cancelled", "skipped", "", "in_progress"):
            passed, calls = self.finish(outcome=result)
            self.assertFalse(passed)
            self.assertTrue(all(call.args[3] == "failure" for call in calls))

    def test_developer_scan_never_changes_early_admission(self):
        for result in ("success", "failure", "cancelled", "skipped"):
            passed, calls = self.finish(trusted=True, outcome=result)
            self.assertEqual(passed, result == "success")
            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0].args[2], security.SCAN)

    def test_scan_finishes_after_merge_without_resetting_admission(self):
        current = pull_request()
        current.update(state="closed", merged=True)
        for trusted in (True, False):
            passed, calls = self.finish(trusted=trusted, current=current)
            self.assertTrue(passed)
            self.assertEqual([call.args[2] for call in calls], [security.SCAN])

    def test_stale_scan_does_not_admit_changed_head_base_or_closed_pr(self):
        for field in ("head", "base", "state"):
            current = pull_request()
            if field == "state":
                current[field] = "closed"
            else:
                current[field]["sha"] = "e" * 40
            passed, calls = self.finish(current=current)
            self.assertFalse(passed)
            self.assertEqual([call.args[2] for call in calls], [security.SCAN])


class AdmissionWaitTests(unittest.TestCase):
    def wait(self, states=("success",), current=None, run=None, target=None, timeout=3000):
        requests = []
        outcomes = iter(states)
        current = current or pull_request()
        event = {"number": 7, "pull_request": {"head": {"sha": HEAD}, "base": {"sha": BASE}}}

        def api(path):
            requests.append(path)
            if "/status?" in path:
                state = next(outcomes)
                return {"statuses": [] if state is None else [{"context": security.ADMISSION, "state": state,
                        "target_url": target or "https://github.com/developer/fds/actions/runs/123"}]}
            if "/actions/runs/" in path:
                return run or {"event": "pull_request_target", "path": ".github/workflows/pr-security.yml"}
            if "/files?" in path:
                return [{"filename": "Source/main.f90"}]
            return current

        with tempfile.TemporaryDirectory() as folder:
            event_path = Path(folder) / "event.json"
            event_path.write_text(json.dumps(event))
            with patch.dict(os.environ, {"GITHUB_EVENT_PATH": str(event_path), "GITHUB_REPOSITORY": "developer/fds",
                                         "GITHUB_SHA": MERGE}), \
                    patch.object(security, "api", side_effect=api), \
                    patch.object(security, "outputs") as output, patch.object(security.time, "sleep") as sleep:
                security.wait_for_admission(timeout=timeout)
                return requests, output, sleep

    def test_waits_until_exact_merge_is_admitted(self):
        requests, output, sleep = self.wait(states=(None, "pending", "success"))
        queried = [path for path in requests if "/status?" in path]
        self.assertEqual(queried, [f"/repos/developer/fds/commits/{MERGE}/status?per_page=100"] * 3)
        self.assertEqual(sleep.call_count, 2)
        output.assert_called_once_with({"conventional": "true", "cmake": "true"})

    def test_failure_error_and_timeout_never_release_builds(self):
        for state in ("failure", "error"):
            with self.subTest(state=state), self.assertRaisesRegex(RuntimeError, "blocked"):
                self.wait(states=(state,))
        with self.assertRaisesRegex(RuntimeError, "Timed out"):
            self.wait(timeout=0)

    def test_status_must_reference_trusted_workflow_in_receiving_repository(self):
        for target in ("https://github.com/outsider/fds/actions/runs/123", "https://example.com/123"):
            with self.assertRaisesRegex(RuntimeError, "trusted workflow"):
                self.wait(target=target)
        for run in ({"event": "pull_request", "path": ".github/workflows/pr-security.yml"},
                    {"event": "pull_request_target", "path": ".github/workflows/other.yml"}):
            with self.assertRaisesRegex(RuntimeError, "security workflow"):
                self.wait(run=run)

    def test_changed_or_closed_pr_never_starts_tests(self):
        for field in ("head", "base", "state"):
            current = pull_request()
            if field == "state":
                current[field] = "closed"
            else:
                current[field]["sha"] = "d" * 40
            with self.subTest(field=field), self.assertRaisesRegex(RuntimeError, "stale revision"):
                self.wait(current=current)


class ArchiveTests(unittest.TestCase):
    def archive(self, folder, files):
        path = folder / "snapshot.tar.gz"
        with tarfile.open(path, "w:gz") as archive:
            for name, content, kind in files:
                info = tarfile.TarInfo(name)
                info.type = kind
                if kind == tarfile.REGTYPE:
                    info.size = len(content)
                    info.mode = 0o777
                    archive.addfile(info, io.BytesIO(content))
                else:
                    info.linkname = "/etc/passwd"
                    archive.addfile(info)
        return path

    def extract(self, files):
        with tempfile.TemporaryDirectory() as name:
            folder = Path(name)
            destination = folder / "source"
            destination.mkdir()
            count = security.extract_snapshot(self.archive(folder, files), destination)
            return count, {str(p.relative_to(destination)): p.read_bytes() for p in destination.rglob("*") if p.is_file()}

    def test_all_extensions_and_hidden_files_are_preserved_without_execution(self):
        files = [("root/.github/workflows/evil.yml", b"test", tarfile.REGTYPE),
                 ("root/run.sh", b"this must never execute", tarfile.REGTYPE),
                 ("root/Manuals/test.pdf", b"pdf test", tarfile.REGTYPE)]
        count, extracted = self.extract(files)
        self.assertEqual(count, 3)
        self.assertEqual(extracted["run.sh"], b"this must never execute")

    def test_unsafe_paths_and_duplicate_entries_fail_closed(self):
        for path in ("../escape", "/root/absolute", "root/../escape", "root/a\n::error::oops", "root/a\\b"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                self.extract([(path, b"data", tarfile.REGTYPE)])
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.extract([("root/same", b"1", tarfile.REGTYPE)] * 2)

    def test_links_devices_submodules_and_lfs_do_not_get_a_clean_result(self):
        for kind in (tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.CHRTYPE, tarfile.FIFOTYPE):
            with self.subTest(kind=kind), self.assertRaises(ValueError):
                self.extract([("root/link", b"", kind)])
        for name, content in (("root/.gitmodules", b"[submodule]"), ("root/large.bin", security.LFS_POINTER + b"\n")):
            with self.subTest(name=name), self.assertRaises(ValueError):
                self.extract([(name, content, tarfile.REGTYPE)])

    def test_oversized_snapshot_fails_instead_of_skipping(self):
        files = [("root/a", b"12345", tarfile.REGTYPE)]
        for limit in ("MAX_FILE", "MAX_TOTAL", "MAX_FILES"):
            with self.subTest(limit=limit), patch.object(security, limit, 0), self.assertRaises(ValueError):
                self.extract(files)

    def test_empty_snapshot_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "Empty"):
            self.extract([])


class ScannerTests(unittest.TestCase):
    def test_only_successful_complete_scans_pass(self):
        cases = [(0, "Scanned files: 1\n", True), (1, "Scanned files: 1\n", False),
                 (2, "Scanned files: 1\n", False), (0, "WARNING: skipped\nScanned files: 1\n", False),
                 (0, "Scanned files: 0\n", False), (0, "", False)]
        with tempfile.TemporaryDirectory() as name:
            path = Path(name)
            (path / "source.f90").write_text("program test\nend program\n")
            for code, output, expected in cases:
                with self.subTest(code=code, output=output), patch.object(security.subprocess, "run", return_value=subprocess.CompletedProcess([], code, output)):
                    if expected:
                        security.clamav(path)
                    else:
                        with self.assertRaises(RuntimeError):
                            security.clamav(path)

    def test_timeout_does_not_pass(self):
        with tempfile.TemporaryDirectory() as name, patch.object(security.subprocess, "run", side_effect=subprocess.TimeoutExpired("clamscan", 1200)):
            with self.assertRaises(subprocess.TimeoutExpired):
                security.clamav(Path(name))

    def test_repository_and_sha_are_validated_before_download(self):
        with patch.object(security.urllib.request, "urlopen") as request:
            for repo, sha in (("../bad", HEAD), ("dev/fds", "master"), ("dev/fds", HEAD + "\n")):
                with self.subTest(repo=repo, sha=sha), self.assertRaises(ValueError):
                    security.download(repo, sha, Path("unused"))
            request.assert_not_called()


class ManifestTests(unittest.TestCase):
    def entry(self, name="source.f90", content=b"source"):
        digest = hashlib.sha1(f"blob {len(content)}\0".encode() + content).hexdigest()
        return {"path": name, "type": "blob", "mode": "100644", "size": len(content), "sha": digest}

    def manifest(self, entries, truncated=False):
        with patch.object(security, "api", side_effect=[
            {"tree": {"sha": BASE}}, {"sha": BASE, "truncated": truncated, "tree": entries},
        ]):
            return security.snapshot_manifest("developer/fds", HEAD)

    def test_git_tree_rejects_truncation_and_hidden_gitlinks(self):
        with self.assertRaises(ValueError):
            self.manifest([self.entry()], truncated=True)
        for kind, mode in (("commit", "160000"), ("blob", "120000")):
            entry = self.entry()
            entry.update(type=kind, mode=mode)
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                self.manifest([entry])

    def test_archive_must_match_every_git_blob(self):
        manifest = self.manifest([self.entry()])
        with tempfile.TemporaryDirectory() as name:
            destination = Path(name)
            path = destination / "source.f90"
            path.write_bytes(b"source")
            security.verify_snapshot(destination, manifest)
            path.write_bytes(b"change")
            with self.assertRaisesRegex(ValueError, "differs"):
                security.verify_snapshot(destination, manifest)
            path.unlink()
            with self.assertRaisesRegex(ValueError, "omits"):
                security.verify_snapshot(destination, manifest)
            path.write_bytes(b"source")
            (destination / "extra").write_bytes(b"extra")
            with self.assertRaisesRegex(ValueError, "adds"):
                security.verify_snapshot(destination, manifest)


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Ruby/Psych is available on GitHub's standard Ubuntu runner and macOS.
        # Parsing, rather than text matching alone, also catches indentation bugs.
        script = ('require "yaml"; require "json"; '
                  'puts JSON.generate(Dir[ARGV[0]].to_h {|f| [File.basename(f), YAML.safe_load(File.read(f))]})')
        result = subprocess.run(["ruby", "-e", script, str(ROOT / ".github/workflows/*.yml")],
                                check=True, capture_output=True, text=True)
        cls.workflows = json.loads(result.stdout)

    def test_builds_use_pull_request_and_security_never_calls_builds(self):
        builds = self.workflows["pr-builds.yml"]
        gate = self.workflows["pr-security.yml"]
        self.assertIn("pull_request", builds.get("on", builds.get("true")))
        self.assertIn("pull_request_target", gate.get("on", gate.get("true")))
        self.assertEqual(set(gate["jobs"]), {"author", "scan", "result"})
        for job in gate["jobs"].values():
            self.assertNotIn("uses", job)
        self.assertNotIn("virus-scan.yml", self.workflows)

    def test_no_pr_code_is_checked_out_before_admission(self):
        gate = self.workflows["pr-security.yml"]
        for job in gate["jobs"].values():
            for step in job["steps"]:
                if step.get("uses", "").startswith("actions/checkout@"):
                    self.assertEqual(step["with"]["ref"], "${{ github.workflow_sha }}")
                    self.assertIs(step["with"]["persist-credentials"], False)
                    self.assertIs(step["with"]["allow-unsafe-checkout"], True)
        wait = self.workflows["pr-builds.yml"]["jobs"]["admission"]["steps"]
        self.assertEqual(wait[0]["with"]["ref"], "${{ github.event.pull_request.base.sha }}")
        self.assertIs(wait[0]["with"]["persist-credentials"], False)
        self.assertIn("pr_security.py wait", wait[1]["run"])

    def test_every_build_depends_on_admission_and_uses_exact_event_merge(self):
        jobs = self.workflows["pr-builds.yml"]["jobs"]
        for name in ("linux", "macos", "windows", "cmake", "line-endings"):
            with self.subTest(job=name):
                self.assertEqual(jobs[name]["needs"], "admission")
                self.assertNotIn("always()", jobs[name].get("if", ""))
                self.assertEqual(jobs[name]["with"]["ref"], "${{ github.sha }}")
                self.assertEqual(jobs[name]["permissions"], {"contents": "read"})
                self.assertNotIn("secrets", jobs[name])

    def test_only_author_and_report_jobs_have_write_permissions(self):
        for name, workflow in self.workflows.items():
            for job_name, job in workflow["jobs"].items():
                permissions = job.get("permissions", workflow.get("permissions", {}))
                writes = [key for key, value in permissions.items() if value == "write"]
                if name == "pr-security.yml" and job_name in ("author", "result"):
                    self.assertEqual(writes, ["statuses"])
                else:
                    self.assertEqual(writes, [])
                if name != "pr-security.yml":
                    for step in job.get("steps", []):
                        self.assertNotIn("allow-unsafe-checkout", step.get("with", {}))

    def test_developer_admission_does_not_depend_on_scan_and_failures_are_reported(self):
        jobs = self.workflows["pr-security.yml"]["jobs"]
        self.assertNotIn("needs", jobs["author"])
        self.assertEqual(jobs["scan"]["needs"], "author")
        self.assertEqual(jobs["scan"]["permissions"], {"contents": "read"})
        self.assertIn("sudo freshclam", jobs["scan"]["steps"][1]["run"])
        self.assertEqual(jobs["result"]["needs"], ["author", "scan"])
        self.assertIn("always()", jobs["result"]["if"])
        self.assertEqual(jobs["result"]["steps"][-1]["env"]["SCAN_RESULT"], "${{ needs.scan.result }}")

    def test_build_summary_fails_for_failed_blocked_or_cancelled_required_jobs(self):
        job = self.workflows["pr-builds.yml"]["jobs"]["result"]
        self.assertEqual(job["name"], "FDS / PR builds")
        self.assertEqual(job["if"], "always()")
        script = job["steps"][0]["run"].split("\n", 1)[1].rsplit("PYTHON", 1)[0]
        results = {name: {"result": "success"} for name in job["needs"]}
        results["admission"]["outputs"] = {"conventional": "true", "cmake": "true"}

        def run(value):
            return subprocess.run([sys.executable, "-c", script],
                                  env={**os.environ, "RESULTS": json.dumps(value)}, capture_output=True).returncode

        self.assertEqual(run(results), 0)
        for name in results:
            for state in ("failure", "cancelled", "skipped"):
                changed = copy.deepcopy(results)
                changed[name]["result"] = state
                with self.subTest(job=name, state=state):
                    self.assertNotEqual(run(changed), 0)
        results["admission"]["outputs"] = {"conventional": "false", "cmake": "false"}
        for name in ("linux", "macos", "windows", "cmake"):
            results[name]["result"] = "skipped"
        self.assertEqual(run(results), 0)


if __name__ == "__main__":
    unittest.main()
