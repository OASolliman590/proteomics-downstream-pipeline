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


SUBPROCESS_CALLS = ("run", "check_output", "Popen")
TEMPFILE_CALLS = ("NamedTemporaryFile", "TemporaryFile", "SpooledTemporaryFile")


def _tempfile_mode(call: ast.Call) -> str:
    for keyword in call.keywords:
        if keyword.arg == "mode" and isinstance(keyword.value, ast.Constant):
            return str(keyword.value.value)
    if call.args and isinstance(call.args[0], ast.Constant) and isinstance(call.args[0].value, str):
        return call.args[0].value
    return "w+b"   # the tempfile default is binary


def violations(source: str, label: str) -> list[str]:
    found = []
    tree = ast.parse(source)
    # review follow-up 2026-10-03: names imported directly (`from subprocess import run as r`, `from tempfile import NamedTemporaryFile`)
    imported: dict[str, tuple[str, str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module in ("subprocess", "tempfile"):
            for alias in node.names:
                imported[alias.asname or alias.name] = (node.module, alias.name)
    for node in ast.walk(tree):
        if isinstance(node, ast.List) and any(isinstance(e, ast.Constant) and e.value == "-e" for e in node.elts) and any(
                isinstance(e, ast.Constant) and e.value == "--vanilla" for e in node.elts):
            found.append(f"{label}:{node.lineno}: `Rscript -e` argv (use runtime.run_r_code)")
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        has_encoding = any(k.arg == "encoding" for k in node.keywords)
        receiver = getattr(getattr(func, "value", None), "id", None)
        module, original = imported.get(name, (receiver, name)) if isinstance(func, ast.Name) else (receiver, name)
        if name in ("read_text", "write_text") and isinstance(func, ast.Attribute) and not has_encoding:
            if not (name == "write_text" and len(node.args) > 1):
                found.append(f"{label}:{node.lineno}: {name}() without encoding")
        elif name == "open" and receiver != "os" and not has_encoding and "b" not in _mode(node):
            found.append(f"{label}:{node.lineno}: text-mode open() without encoding")
        elif original in SUBPROCESS_CALLS and module == "subprocess" and not has_encoding and any(
                k.arg in ("text", "universal_newlines") and isinstance(k.value, ast.Constant) and k.value.value for k in node.keywords):
            found.append(f"{label}:{node.lineno}: text-mode subprocess.{original}() without encoding")
        elif original in TEMPFILE_CALLS and module == "tempfile" and not has_encoding and "b" not in _tempfile_mode(node):
            found.append(f"{label}:{node.lineno}: text-mode tempfile.{original}() without encoding")
    return found


# R sources (review follow-up 2026-10-03): reads declare UTF-8 explicitly instead of relying on R >= 4.2's native UTF-8 on
# Windows, and text artifacts are written through the LF writers (.pc_write_tsv/.pc_write_json/.pc_write_lf), because a
# text-mode connection writes CRLF on Windows and changes artifact bytes.
R_READS = ("read.delim", "read.table", "read.csv", "readLines", "scan")


def _r_calls(source: str, name: str):
    import re
    for match in re.finditer(r"(?<![A-Za-z0-9_.])" + re.escape(name) + r"\s*\(", source):
        depth, i = 0, match.end() - 1
        while i < len(source):
            depth += {"(": 1, ")": -1}.get(source[i], 0)
            if depth == 0:
                break
            i += 1
        yield source.count("\n", 0, match.start()) + 1, source[match.start():i + 1]


def r_violations(source: str, label: str) -> list[str]:
    found = []
    for name in R_READS:
        for line, call in _r_calls(source, name):
            if "encoding" not in call.lower():
                found.append(f"{label}:{line}: {name}() without an explicit UTF-8 encoding")
    for line, call in _r_calls(source, "jsonlite::write_json"):
        found.append(f"{label}:{line}: jsonlite::write_json() (text mode; use .pc_write_json)")
    for line, call in _r_calls(source, "write.table"):
        if 'eol = "\\n"' not in call and 'eol="\\n"' not in call:
            found.append(f"{label}:{line}: write.table() without eol = \"\\n\"")
    return found


def test_every_r_text_io_names_utf8_and_lf():
    problems = []
    for path in sorted((ROOT / "r" / "proteomicsCore" / "R").glob("*.R")):
        problems += r_violations(path.read_text(encoding="utf-8"), path.relative_to(ROOT).as_posix())
    assert problems == [], "\n".join(problems)


def test_r_lint_detects_locale_dependent_io():
    bad = 'x <- read.delim(path, sep = "\\t")\njsonlite::write_json(v, path)\nwrite.table(d, file = path, sep = "\\t")\n'
    good = 'x <- read.delim(path, sep = "\\t", encoding = "UTF-8")\n.pc_write_json(v, path)\nwrite.table(d, con, sep = "\\t", eol = "\\n")\n'
    assert len(r_violations(bad, "bad")) == 3 and r_violations(good, "good") == []


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
    # review follow-up 2026-10-03: direct imports, tempfile text modes and `Rscript -e` argv lists
    bad = ("from subprocess import run, check_output as co\nimport tempfile\nfrom tempfile import NamedTemporaryFile as N\n"
           "run(['x'], text=True)\nco(['x'], universal_newlines=True)\ntempfile.NamedTemporaryFile('w')\nN(mode='w+')\ntempfile.TemporaryFile(mode='r+')\n"
           "subprocess.run([r, '--vanilla', '-e', 'quit()'])\n")
    good = ("from subprocess import run\nimport tempfile\nfrom tempfile import NamedTemporaryFile as N\nrun(['x'], text=True, encoding='utf-8')\n"
            "tempfile.NamedTemporaryFile()\nN(mode='w', encoding='utf-8')\ntempfile.TemporaryFile('w+b')\nrun([r, '--vanilla', script])\n")
    assert len(violations(bad, "bad")) == 6 and violations(good, "good") == []


# Review 2026-10-05 (minor 6): files that write hashed release/regression artifacts must write LF bytes on every
# platform (D-42). A text-mode write without newline="\n" translates "\n" to os.linesep (CRLF on Windows).
NEWLINE_STRICT = [ROOT / "scripts" / "maintained" / "build_release.py", ROOT / "src" / "proteomics_pipeline" / "legacy_service.py"]


def newline_violations(source: str, label: str) -> list[str]:
    found = []
    for node in ast.walk(ast.parse(source)):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", None)
        has_newline = any(k.arg == "newline" for k in node.keywords)
        if name == "write_text" and isinstance(func, ast.Attribute) and not has_newline:
            found.append(f"{label}:{node.lineno}: write_text() without newline=\"\\n\" (CRLF on Windows)")
        elif name == "open" and any(c in _mode(node) for c in "wax") and "b" not in _mode(node) and not has_newline:
            found.append(f"{label}:{node.lineno}: text-mode open() for writing without newline=\"\\n\"")
    return found


def test_release_artifact_writers_use_lf_on_every_platform():
    problems = []
    for path in NEWLINE_STRICT:
        problems += newline_violations(path.read_text(encoding="utf-8"), path.relative_to(ROOT).as_posix())
    assert problems == [], "\n".join(problems)


def test_newline_lint_detects_platform_dependent_writes():
    bad = "from pathlib import Path\nPath('m').write_text('x\\n', encoding='utf-8')\nopen('r', 'w', encoding='utf-8')\n"
    good = "from pathlib import Path\nPath('m').write_bytes(b'x\\n')\nPath('m').write_text('x\\n', encoding='utf-8', newline='\\n')\nopen('r', 'w', encoding='utf-8', newline='\\n')\nopen('r', encoding='utf-8')\n"
    assert len(newline_violations(bad, "bad")) == 2 and newline_violations(good, "good") == []
