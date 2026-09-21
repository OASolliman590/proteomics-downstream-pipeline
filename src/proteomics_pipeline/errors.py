class ProteomicsError(Exception):
    def __init__(self, code: str, message: str, pointer: str | None = None, *, exit_code: int = 2):
        self.code, self.message, self.pointer, self.exit_code = code, message, pointer, exit_code
        super().__init__(message)
    def as_dict(self) -> dict:
        result = {"code": self.code, "message": self.message}
        if self.pointer is not None: result["pointer"] = self.pointer
        return result
class ConfigurationError(ProteomicsError): pass
class CapabilityError(ProteomicsError):
    def __init__(self, code: str, message: str): super().__init__(code, message, exit_code=3)
class IntegrityError(ProteomicsError):
    def __init__(self, message: str, code: str = "E_INTEGRITY"): super().__init__(code, message, exit_code=5)
class CollisionError(IntegrityError):
    def __init__(self, message: str = "an exclusive run writer already holds the lock"): super().__init__(message, code="E_RUN_COLLISION")
