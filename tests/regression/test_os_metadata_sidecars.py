"""Operating-system metadata sidecars are never pipeline files (reconciliation 2026-10-18, item 4).

macOS writes an AppleDouble `._<name>` file next to every file on exFAT/FAT/SMB volumes (and `.DS_Store` in folders);
Windows Explorer writes `Thumbs.db` and `desktop.ini`. A sidecar named `._figure.png` or `._module.py` must not be
registered as a figure artifact, hashed into a manifest or code identity, or parsed as a source table.
"""
from __future__ import annotations

from pathlib import Path

from proteomics_pipeline.figures.run_catalogue import discover
from proteomics_pipeline.provenance import file_manifest, is_os_metadata
from proteomics_pipeline.workflow import code_manifest

APPLEDOUBLE = b"\x00\x05\x16\x07\x00\x02\x00\x00Mac OS X        \x00\x02\xff\xfe\x80\x81"   # not valid UTF-8


def _sidecars(folder: Path, name: str) -> None:
    (folder / f"._{name}").write_bytes(APPLEDOUBLE)
    (folder / ".DS_Store").write_bytes(b"\x00\x00\x00\x01Bud1\xff")


def test_sidecar_names_are_recognised_and_real_names_are_not():
    for name in ("._figure.png", "._module.py", ".DS_Store", "Thumbs.db", "desktop.ini"):
        assert is_os_metadata(Path("run") / name), name
    for name in ("figure.png", "_private.py", ".stage-01", "a._b.tsv", "DS_Store"):
        assert not is_os_metadata(Path("run") / name), name


def test_file_manifest_ignores_sidecars(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "table.tsv").write_text("a\tb\n1\t2\n", encoding="utf-8")
    before = file_manifest(tmp_path)
    _sidecars(tmp_path / "sub", "table.tsv")
    (tmp_path / "Thumbs.db").write_bytes(b"\xd0\xcf\x11\xe0")
    assert file_manifest(tmp_path) == before == [{"relative_path": "sub/table.tsv", "sha256": before[0]["sha256"]}]


def test_code_identity_is_unchanged_by_sidecars(tmp_path):
    package = tmp_path / "proteomics_pipeline"
    package.mkdir()
    (package / "module.py").write_text("X = 1\n", encoding="utf-8")
    before = code_manifest(package, tmp_path)
    _sidecars(package, "module.py")
    assert code_manifest(package, tmp_path) == before


def test_run_catalogue_never_reads_a_sidecar_as_a_source_table(tmp_path):
    sources = tmp_path / "stages" / "differential" / "figure_sources"
    sources.mkdir(parents=True)
    _sidecars(sources, "volcano_primary.tsv")
    found, skipped = discover(tmp_path)
    assert found == [] and skipped == []
