"""Optional web layer: account-scoped genome upload + a per-user traits API.

The FastAPI stack is an optional dependency (``pip install -e '.[web]'``); import
:func:`create_app` lazily so the rest of the package works without it.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

__all__ = ["create_app"]

if TYPE_CHECKING:
    from .app import create_app


def create_app(*args, **kwargs):
    """Build the FastAPI app (imports the web stack on first call)."""
    try:
        from .app import create_app as _create_app
    except ImportError as exc:  # pragma: no cover - exercised via the CLI
        raise RuntimeError(
            "The web layer needs the optional web dependencies: "
            "pip install -e '.[web]'"
        ) from exc
    return _create_app(*args, **kwargs)
