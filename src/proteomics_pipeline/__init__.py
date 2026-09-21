"""Maintained orchestration foundation for the proteomics pipeline."""
__version__ = "0.1.0.dev0"

def capabilities():
    """Advertise the one foundation handler owned by R01."""
    return [{"id": "foundation.io_roundtrip", "dependencies": [],
             "required_r_packages": ["jsonlite", "openssl", "proteomicsCore"]}]

def io_roundtrip(request):
    """Lazy bridge to the orchestrator's foundation execution seam."""
    from .runtime import execute_stage
    return execute_stage(request)

__all__ = ["__version__", "capabilities", "io_roundtrip"]
