"""Context-managed webcam capture (OpenCV).

Kept deliberately thin: it only yields frames.  Face detection lives in
:mod:`magic_mirror.detect`, and decoration in :mod:`magic_mirror.render`.
"""

from __future__ import annotations

import numpy as np

try:  # pragma: no cover - import guard exercised only by environment
    import cv2
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "The live mirror needs OpenCV. Install it with "
        "`pip install 'genetic-reflexions[camera]'`, or use the poster/list "
        "commands which have no OpenCV dependency."
    ) from exc


class Camera:
    """Context-managed webcam that yields mirror-flipped RGB frames."""

    def __init__(self, device_id: int, mirror: bool = True) -> None:
        self.device_id = device_id
        self.mirror = mirror
        self._capture: cv2.VideoCapture | None = None

    def __enter__(self) -> Camera:
        self._capture = cv2.VideoCapture(self.device_id)
        if not self._capture.isOpened():
            raise RuntimeError(
                f"Could not open camera device {self.device_id}. "
                "Check it is connected and not in use by another app."
            )
        return self

    def read_rgb(self) -> np.ndarray | None:
        """Return the next frame as RGB, mirror-flipped, or ``None`` on failure."""
        assert self._capture is not None
        ok, frame_bgr = self._capture.read()
        if not ok or frame_bgr is None:
            return None
        if self.mirror:
            frame_bgr = cv2.flip(frame_bgr, 1)
        return cv2.cvtColor(frame_bgr, cv2.COLOR_BGR2RGB)

    def __exit__(self, *_exc) -> None:
        if self._capture is not None:
            self._capture.release()
