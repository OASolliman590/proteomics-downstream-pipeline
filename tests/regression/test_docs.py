"""R12 V115: usage and methods documentation agrees with the implemented handlers (FR-115, SM25).

Oracles: the runtime's own capability discovery, the configuration schema enums, and the reason codes found in the
maintained sources (an independent regex scan written here). The documented commands themselves are executed by
R10b V100 (tests/integration/test_full_report.py) on the synthetic examples.
"""
from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

from proteomics_pipeline.runtime import CAPABILITY_MAP, capabilities as discovered_capabilities

ROOT = Path(__file__).resolve().parents[2]
METHODS = ROOT / "docs" / "methods" / "README.md"
CODES = ROOT / "docs" / "methods" / "reason-codes.md"
RELEASE = ROOT / "docs" / "user-guide" / "release.md"
PUBLIC_DOCS = [ROOT / "README.md", ROOT / "CHANGELOG.md", ROOT / "docs" / "user-guide" / "usage.md", METHODS, RELEASE]
spec = importlib.util.spec_from_file_location("build_release", ROOT / "scripts" / "maintained" / "build_release.py")
BR = importlib.util.module_from_spec(spec); spec.loader.exec_module(BR)
# D-63: the Maintainer selected PolyForm Noncommercial 1.0.0 (docs/validation/013-release/operator-decisions.json).
# Any claim of a different license is unsupported.
OTHER_LICENSE_CLAIM = re.compile(r"(?i)\bMIT License\b|Apache License|GNU (?:Affero |Lesser )?General Public License|BSD [0-9]-Clause|CC-BY|Creative Commons|SPDX-License-Identifier")
RECORDED_LICENSE = "PolyForm Noncommercial License 1.0.0"
REQUIRED_NOTICE = "Required Notice: Copyright (c) 2026 Omar A. Solliman"


def _documented_capabilities() -> set[str]:
    table = METHODS.read_text(encoding="utf-8").split("## Capabilities", 1)[1].split("\n## ", 1)[0]
    return {c for row in table.splitlines() if row.startswith("| `") for c in re.findall(r"`([a-z_.]+)`", row.split("|")[1])}


def test_v115_every_implemented_capability_is_documented_and_nothing_else():
    implemented = {c["id"] for c in discovered_capabilities() if c["implemented"]}
    documented = _documented_capabilities()
    assert documented == implemented, (documented ^ implemented)
    assert documented <= set(CAPABILITY_MAP)


def test_v115_documented_declarations_equal_the_schema_enums():
    schema = json.loads((ROOT / "src" / "proteomics_pipeline" / "schemas" / "analysis.schema.json").read_text(encoding="utf-8"))
    text = METHODS.read_text(encoding="utf-8")
    enums = [schema["properties"]["assay"]["enum"], schema["properties"]["input"]["properties"]["format"]["enum"], schema["properties"]["source_scale"]["enum"],
             schema["$defs"]["model"]["properties"]["engine"]["enum"], schema["$defs"]["blocking"]["properties"]["mode"]["enum"],
             schema["$defs"]["coverage"]["properties"]["policy"]["enum"], schema["$defs"]["sensitivity"]["properties"]["method"]["enum"],
             schema["properties"]["preprocessing"]["properties"]["normalization"]["enum"],
             schema["$defs"]["post_de"]["properties"]["biomarker"]["properties"]["classifier"]["enum"],
             schema["$defs"]["post_de"]["properties"]["biomarker"]["properties"]["selection"]["properties"]["method"]["enum"],
             schema["$defs"]["post_de"]["properties"]["networks"]["properties"]["coabundance"]["properties"]["rule"]["enum"]]
    missing = [v for enum in enums for v in enum if f"`{v}`" not in text]
    assert missing == []


def test_v115_reason_code_catalogue_matches_the_sources_both_ways():
    scanned = set()
    for path in list((ROOT / "src" / "proteomics_pipeline").rglob("*.py")) + list((ROOT / "r" / "proteomicsCore" / "R").glob("*.R")):
        scanned |= set(re.findall(r"\b([EW]_[A-Z][A-Z0-9_]+)\b", path.read_text(encoding="utf-8")))
    documented = set(re.findall(r"^\| `([EW]_[A-Z0-9_]+)` \|", CODES.read_text(encoding="utf-8"), re.M))
    assert documented == scanned, (sorted(documented - scanned), sorted(scanned - documented))
    assert CODES.read_text(encoding="utf-8") == BR.codes_markdown(ROOT)          # regenerated, not hand-edited


def test_v115_negative_license_claims_match_the_recorded_decision_and_no_v1_completion_claim():
    decision = json.loads((ROOT / BR.DECISION_RECORD).read_text(encoding="utf-8"))
    assert decision["license"] == "PolyForm-Noncommercial-1.0.0" and decision["license_file"] == "LICENSE.md"
    for path in PUBLIC_DOCS:
        assert not OTHER_LICENSE_CLAIM.search(path.read_text(encoding="utf-8")), path
    for path in (ROOT / "README.md", METHODS, RELEASE, ROOT / "CHANGELOG.md"):
        text = path.read_text(encoding="utf-8")
        assert RECORDED_LICENSE in text and "No license has been selected" not in text, path
    for path in (METHODS, RELEASE):
        assert "not a v1.0 release" in path.read_text(encoding="utf-8"), path
    license_files = sorted(p.name for p in ROOT.iterdir() if p.name.upper().startswith(("LICENSE", "LICENCE", "COPYING")))
    assert license_files == ["LICENSE.md"]
    text = (ROOT / "LICENSE.md").read_text(encoding="utf-8")
    assert text.startswith("# " + RECORDED_LICENSE) and REQUIRED_NOTICE in text
    assert (ROOT / "r" / "proteomicsCore" / "LICENSE").read_text(encoding="utf-8") == text
    assert "License: file LICENSE" in (ROOT / "r" / "proteomicsCore" / "DESCRIPTION").read_text(encoding="utf-8")


def test_v115_negative_detector_catches_an_unsupported_claim(tmp_path):
    assert OTHER_LICENSE_CLAIM.search("This project is licensed under the MIT License.")
    assert not OTHER_LICENSE_CLAIM.search("Licensed under the " + RECORDED_LICENSE + ".")
    fake = "| `deep_learning_engine` | R99 | not implemented |"
    table = METHODS.read_text(encoding="utf-8").replace("## Capabilities\n", "## Capabilities\n" + fake + "\n", 1)
    rows = {c for row in table.split("## Capabilities", 1)[1].split("\n## ", 1)[0].splitlines() if row.startswith("| `") for c in re.findall(r"`([a-z_.]+)`", row.split("|")[1])}
    assert "deep_learning_engine" in rows and "deep_learning_engine" not in {c["id"] for c in discovered_capabilities() if c["implemented"]}
