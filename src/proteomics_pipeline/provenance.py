from __future__ import annotations
import hashlib, json
from pathlib import Path
from typing import Any
def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
_PLAN_EXCLUDED_KEYS=frozenset({"plan_hash","started_at","finished_at","created_at","generated_at","observed_at","timestamp"})
def _plan_value(value: Any) -> Any:
    if isinstance(value,dict): return {key:_plan_value(item) for key,item in value.items() if key not in _PLAN_EXCLUDED_KEYS}
    if isinstance(value,list): return [_plan_value(item) for item in value]
    return value
def canonical_json_sha256(value: Any) -> str:
    """Hash the frozen plan representation (excluding mutable observations)."""
    return hashlib.sha256(canonical_plan_bytes(value)).hexdigest()
def canonical_config_sha256(value: Any) -> str:
    """Hash a raw schema-validated configuration without plan exclusions."""
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()
def canonical_plan_bytes(value: dict[str,Any]) -> bytes: return canonical_json_bytes(_plan_value(value))
def canonical_plan_sha256(value: dict[str,Any]) -> str: return canonical_json_sha256(value)
def sha256_bytes(value: bytes) -> str: return hashlib.sha256(value).hexdigest()
def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""): digest.update(chunk)
    return digest.hexdigest()
# Decision D-42 (2026-10-03, review follow-up): platform-independent content identity.  A Windows checkout
# (core.autocrlf / text=auto) or a Windows R writer turns LF into CRLF without changing the content, which made the
# plan hash differ between CI platforms (run 37087074374).  Text content is therefore identified by the SHA-256 of
# its bytes with every CRLF pair replaced by LF; files containing NUL bytes are treated as binary and hashed raw.
CONTENT_HASH_RULE = "sha256 over bytes with CRLF normalised to LF; files containing NUL bytes are hashed raw"
def normalized_content(data: bytes) -> bytes:
    return data if b"\x00" in data else data.replace(b"\r\n", b"\n")
def content_sha256(path: str | Path) -> str:
    """Platform-independent content hash (see CONTENT_HASH_RULE)."""
    return hashlib.sha256(normalized_content(Path(path).read_bytes())).hexdigest()
def file_manifest(root: str | Path) -> list[dict[str, str]]:
    root = Path(root).resolve()
    return [{"relative_path": p.relative_to(root).as_posix(), "sha256": sha256_file(p)} for p in sorted(x for x in root.rglob("*") if x.is_file())]
