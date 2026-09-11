"""Offline regression tests; no real malware/EICAR files or network requests."""

import copy
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
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
    def test_public_organization_members_only(self):
        with patch.object(security, "api", return_value=None) as api:
            self.assertTrue(security.organization_member({"login": "member", "type": "User"}))
            api.assert_called_once_with("/orgs/firemodels/public_members/member", authenticated=False)

    def test_private_unknown_or_unavailable_membership_requires_scan(self):
        for error in (urllib.error.HTTPError("url", 404, "Not Found", {}, None),
                      urllib.error.HTTPError("url", 403, "Rate limited", {}, None),
                      urllib.error.URLError("offline"), TimeoutError(), ValueError()):
            with self.subTest(error=error), patch.object(security, "api", side_effect=error):
                self.assertFalse(security.organization_member({"login": "author", "type": "User"}))

    def test_bots_and_invalid_names_are_not_members(self):
        with patch.object(security, "api") as api:
            for user in ({"login": "dependabot[bot]", "type": "Bot"},
                         {"login": "../member", "type": "User"}, {"login": "member"}):
                self.assertFalse(security.organization_member(user))
            api.assert_not_called()


class PreparationTests(unittest.TestCase):
    def prepare(self, pull=None, memberships=(False, False), event_sha=HEAD):
        pull = pull or pull_request()
        event = {"number": 7, "pull_request": {"head": {"sha": event_sha}},
                 "sender": {"login": "sender", "type": "User"}}

        def api(path, **kwargs):
            if path.endswith("/git/commits/" + MERGE):
                return {"parents": [{"sha": BASE}, {"sha": HEAD}]}
            if "/files?" in path:
                return [{"filename": "Source/main.f90"}]
            return copy.deepcopy(pull)

        with tempfile.TemporaryDirectory() as folder:
            event_path = Path(folder) / "event.json"
            event_path.write_text(json.dumps(event))
            with patch.dict(os.environ, {"GITHUB_EVENT_PATH": str(event_path), "GITHUB_REPOSITORY": "developer/fds"}), \
                    patch.object(security, "api", side_effect=api), \
                    patch.object(security, "status") as status, \
                    patch.object(security, "organization_member", side_effect=memberships), \
                    patch.object(security, "outputs") as output, patch.object(security.time, "sleep"):
                security.prepare()
                self.assertEqual(status.call_args_list[0].args[:2], ("developer/fds", event_sha))
                return output.call_args.args[0]

    def test_fork_owner_or_existing_collaborator_is_not_automatically_trusted(self):
        pull = pull_request()
        pull["author_association"] = "OWNER"
        result = self.prepare(pull)
        self.assertEqual(result["trusted"], "false")
        self.assertEqual(result["merge_sha"], MERGE)

    def test_author_and_sender_must_both_be_organization_members(self):
        self.assertEqual(self.prepare(memberships=(True, True))["trusted"], "true")
        self.assertEqual(self.prepare(memberships=(True, False))["trusted"], "false")
        self.assertEqual(self.prepare(memberships=(False, True))["trusted"], "false")

    def test_stale_event_does_not_admit_a_new_revision(self):
        with self.assertRaisesRegex(ValueError, "stale"):
            self.prepare(event_sha="d" * 40)

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
    def check(self, trusted=False, external="success", member="success", merge=MERGE, current=None, kind="admission"):
        environment = {"GITHUB_REPOSITORY": "developer/fds", "HEAD_SHA": HEAD, "BASE_SHA": BASE,
                       "PR_NUMBER": "7", "MERGE_SHA": merge, "TRUSTED": str(trusted).lower(),
                       "EXTERNAL_RESULT": external, "TRUSTED_RESULT": member}
        with patch.dict(os.environ, environment), patch.object(security, "api", return_value=current or pull_request()), \
                patch.object(security, "status") as status:
            try:
                security.report(kind)
                success = True
            except RuntimeError:
                success = False
            return success, status

    def test_external_scan_must_finish_successfully(self):
        self.assertTrue(self.check(external="success")[0])
        for outcome in ("failure", "cancelled", "skipped", "", "in_progress"):
            with self.subTest(outcome=outcome):
                success, status = self.check(external=outcome)
                self.assertFalse(success)
                self.assertEqual(status.call_args.args[3], "failure")

    def test_member_admission_is_independent_of_background_scan(self):
        for outcome in ("success", "failure", "cancelled", "in_progress"):
            self.assertTrue(self.check(trusted=True, external="skipped", member=outcome)[0])

    def test_member_scan_failure_is_still_reported(self):
        success, status = self.check(trusted=True, external="skipped", member="failure", kind="report")
        self.assertFalse(success)
        self.assertEqual(status.call_args.args[2:4], (security.SCAN, "failure"))

    def test_missing_merge_snapshot_never_admits(self):
        for trusted in (True, False):
            self.assertFalse(self.check(trusted=trusted, merge="")[0])

    def test_changed_head_base_or_closed_pr_does_not_overwrite_new_status(self):
        for field in ("head", "base", "state"):
            current = pull_request()
            if field == "state":
                current[field] = "closed"
            else:
                current[field]["sha"] = "e" * 40
            success, status = self.check(current=current)
            self.assertFalse(success)
            status.assert_not_called()

    def test_build_status_requires_all_selected_checks(self):
        results = {name: {"result": "success"} for name in
                   ("prepare", "admission", "linux", "macos", "windows", "cmake", "line-endings")}
        environment = {"BUILD_RESULTS": json.dumps(results), "CONVENTIONAL": "true", "CMAKE": "true"}
        with patch.dict(os.environ, environment):
            passed, status = self.check(kind="builds")
            self.assertTrue(passed)
            self.assertEqual(status.call_args.args[2], security.BUILDS)
        for name in results:
            changed = copy.deepcopy(results)
            changed[name]["result"] = "skipped"
            with self.subTest(job=name), patch.dict(os.environ, {**environment, "BUILD_RESULTS": json.dumps(changed)}):
                self.assertFalse(self.check(kind="builds")[0])
        for name in ("linux", "macos", "windows", "cmake"):
            results[name]["result"] = "skipped"
        with patch.dict(os.environ, {"BUILD_RESULTS": json.dumps(results), "CONVENTIONAL": "false", "CMAKE": "false"}):
            self.assertTrue(self.check(kind="builds")[0])


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

    def test_pr_code_cannot_directly_trigger_the_existing_workflows(self):
        for name, workflow in self.workflows.items():
            events = workflow.get("on", workflow.get("true"))
            with self.subTest(workflow=name):
                self.assertNotIn("pull_request", events)
                if name in ("pr-security.yml", "virus-scan.yml"):
                    continue
                self.assertIn("workflow_call", events)
                self.assertIn("push", events)
                for job in workflow["jobs"].values():
                    for step in job.get("steps", []):
                        if step.get("uses", "").startswith("actions/checkout@"):
                            self.assertEqual(step["with"]["ref"], "${{ inputs.ref || github.sha }}")
                            self.assertIs(step["with"]["persist-credentials"], False)

    def test_every_pr_build_depends_on_successful_admission(self):
        jobs = self.workflows["pr-security.yml"]["jobs"]
        for name in ("linux", "macos", "windows", "cmake", "line-endings"):
            with self.subTest(job=name):
                self.assertIn("admission", jobs[name]["needs"])
                self.assertIn("needs.admission.result == 'success'", jobs[name]["if"])
                self.assertEqual(jobs[name]["with"]["ref"], "${{ needs.prepare.outputs.merge_sha }}")
                self.assertEqual(jobs[name]["permissions"], {"contents": "read"})
                self.assertNotIn("secrets", jobs[name])

    def test_reusable_workflows_do_not_share_concurrency_groups(self):
        groups = []
        for name in ("linux.yml", "osx.yml", "windows.yml", "cmake.yml"):
            group = self.workflows[name]["concurrency"]["group"]
            self.assertIn("github.event.pull_request.number", group)
            self.assertNotIn("github.workflow", group)
            groups.append(group)
        self.assertEqual(len(set(groups)), len(groups))

    def test_member_scan_is_not_an_admission_dependency(self):
        jobs = self.workflows["pr-security.yml"]["jobs"]
        self.assertEqual(jobs["admission"]["needs"], ["prepare", "scan-external"])
        self.assertIn("scan-member", jobs["scan-result"]["needs"])
        self.assertEqual(jobs["scan-external"]["if"], "needs.prepare.outputs.trusted != 'true'")
        self.assertEqual(jobs["scan-member"]["if"], "needs.prepare.outputs.trusted == 'true'")
        self.assertNotIn("scan-member", jobs["build-result"]["needs"])
        self.assertNotIn("scan-external", jobs["build-result"]["needs"])

    def test_scan_and_identity_checkouts_use_only_trusted_policy(self):
        for name in ("pr-security.yml", "virus-scan.yml"):
            for job in self.workflows[name]["jobs"].values():
                for step in job.get("steps", []):
                    if step.get("uses", "").startswith("actions/checkout@"):
                        self.assertRegex(step["uses"], r"@[0-9a-f]{40}$")
                        self.assertIn(step["with"]["ref"], ("${{ github.workflow_sha }}", "${{ inputs.policy_sha }}"))
                        self.assertIs(step["with"]["persist-credentials"], False)
        scanner = self.workflows["virus-scan.yml"]
        self.assertEqual(scanner["permissions"], {"contents": "read"})
        steps = scanner["jobs"]["clamav"]["steps"]
        self.assertIn("sudo freshclam", steps[1]["run"])
        self.assertEqual(steps[2]["env"]["GH_TOKEN"], "${{ github.token }}")


if __name__ == "__main__":
    unittest.main()
