"""Tests for the read-only external-data archive audit helper."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import tarfile
import zipfile


SCRIPT = (
    Path(__file__).parents[2] / "scripts" / "audit_exoctk_data_archives.py"
)
SPEC = spec_from_file_location("audit_exoctk_data_archives", SCRIPT)
audit = module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_discover_archive_urls_expands_patch_version(tmp_path):
    package = tmp_path / "exoctk"
    package.mkdir()
    (package / "_version.py").write_text("__version__ = '2026.7.1'\n")
    (tmp_path / "urls.py").write_text(
        "URL = 'https://" + "example.test/data{PATCHVER}.tar.gz'\n"
    )

    urls = audit.discover_archive_urls(tmp_path)

    expected_url = "https://" + "example.test/datav2026.7.tar.gz"
    assert urls == {expected_url: ["urls.py"]}


def test_zip_inspection_reports_duplicates_metadata_and_unsafe_paths(tmp_path):
    archive = tmp_path / "synthetic.zip"
    with zipfile.ZipFile(archive, "w") as stream:
        stream.writestr("data/value.txt", "first")
        stream.writestr("data/value.txt", "second")
        stream.writestr("__MACOSX/._value.txt", "metadata")
        stream.writestr("../escape.txt", "unsafe")

    report = audit.inspect_archive(archive)

    assert report["archive_format"] == "zip"
    assert report["member_count"] == 4
    assert report["duplicate_member_paths"] == ["data/value.txt"]
    assert report["metadata_artifacts"] == ["__MACOSX/._value.txt"]
    assert report["unsafe_member_paths"] == ["../escape.txt"]


def test_tar_inspection_does_not_extract_members(tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("audit only\n")
    archive = tmp_path / "synthetic.tar.gz"
    with tarfile.open(archive, "w:gz") as stream:
        stream.add(source, arcname="family/source.txt")
    source.unlink()

    report = audit.inspect_archive(archive)

    assert report["archive_format"] == "tar"
    assert report["top_level_paths"] == {"family": 1}
    assert not source.exists()
