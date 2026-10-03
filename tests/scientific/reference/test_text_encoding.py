"""Portability lint (R11 V107 follow-up, CI run 36982402784).

Windows runners decoded UTF-8 files and R output with the cp1252 code page
(`UnicodeDecodeError: 'charmap' codec ...`, mojibake 'Î©').  Every text-mode
open()/read_text()/write_text() and every text-mode subprocess call in the
maintained code and tests must name encoding="utf-8" explicitly.  The
recovered study runner (scripts/run_pipeline.py) is baseline code and is not
edited.
"""
from __future__ import annotations

import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCOPES = [ROOT / "src", ROOT / "tests", ROOT / "scripts" / "maintained"]


def _mode(call: ast.Call) -> str:
    for keyword in call.keywords:
        if keyword.arg == "mode" and isinstance(keyword.value, ast.Constant):
            return str(keyword.value.value)
    position = 1 if isinstance(call.func, ast.Name) else 0
    if len(call.args) > position and isinstance(call.args[position], ast.Constant) and isinstance(call.args[position].value, str):
        return call.args[position].value
    return "r"


def violations(source: str, label: str) -> list[str]:
    found = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        has_encoding = any(k.arg == "encoding" for k in node.keywords)
        receiver = getattr(getattr(func, "value", None), "id", None)
        if name in ("read_text", "write_text") and isinstance(func, ast.Attribute) and not has_encoding:
            if not (name == "write_text" and len(node.args) > 1):
                found.append(f"{label}:{node.lineno}: {name}() without encoding")
        elif name == "open" and receiver != "os" and not has_encoding and "b" not in _mode(node):
            found.append(f"{label}:{node.lineno}: text-mode open() without encoding")
        elif name in ("run", "check_output", "Popen") and receiver == "subprocess" and not has_encoding and any(
                k.arg in ("text", "universal_newlines") and isinstance(k.value, ast.Constant) and k.value.value for k in node.keywords):
            found.append(f"{label}:{node.lineno}: text-mode subprocess.{name}() without encoding")
    return found


def test_every_text_io_names_utf8():
    problems = []
    for scope in SCOPES:
        for path in sorted(scope.rglob("*.py")):
            problems += violations(path.read_text(encoding="utf-8"), path.relative_to(ROOT).as_posix())
    assert problems == [], "\n".join(problems)


def test_lint_detects_locale_dependent_io():
    bad = "from pathlib import Path\nimport subprocess\nPath('a').read_text()\nopen('b')\nsubprocess.run(['x'], text=True)\n"
    good = "from pathlib import Path\nimport os, subprocess\nPath('a').read_text(encoding='utf-8')\nopen('b', 'rb')\nos.open('c', 0)\nsubprocess.run(['x'], text=True, encoding='utf-8')\n"
    assert len(violations(bad, "bad")) == 3 and violations(good, "good") == []
