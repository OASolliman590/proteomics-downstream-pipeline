from __future__ import annotations
import os, tempfile
from contextlib import contextmanager
from pathlib import Path
def resolve_within(root: str | Path, relative: str | Path) -> Path:
    root_path = Path(root).resolve(); candidate = (root_path / relative).resolve()
    if candidate != root_path and root_path not in candidate.parents: raise ValueError(f"path escapes root: {relative}")
    return candidate
@contextmanager
def temporary_directory(root: str | Path, prefix: str = ".stage-"):
    root_path = Path(root); root_path.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=prefix, dir=root_path) as name: yield Path(name)
def atomic_write_bytes(destination: str | Path, data: bytes) -> None:
    destination = Path(destination); destination.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data); handle.flush(); os.fsync(handle.fileno())
        os.replace(name, destination)
    except Exception:
        try: os.unlink(name)
        except FileNotFoundError: pass
        raise
