"""Python interface for the C++ solar dynamo model."""

try:
    from ._core import Model, Parameters, Snapshot
except ImportError as exc:
    raise ImportError(
        "The exact legacy solver requires the compiled C++ extension. "
        "Install with `pip install -e .`; no approximate fallback is used."
    ) from exc
BACKEND = "C++ exact-v3"

__all__ = ["BACKEND", "Model", "Parameters", "Snapshot"]
__version__ = "0.2.0"
