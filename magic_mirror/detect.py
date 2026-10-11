"""Face detection, newest-first with graceful fallback.

- :class:`YuNetDetector` — OpenCV's modern YuNet DNN detector (accurate, fast).
  The ~230 KB ONNX model is downloaded once and cached on disk.
- :class:`HaarDetector` — the classic Haar cascade, used if YuNet/the model is
  unavailable so the mirror still runs anywhere.
- :class:`SmoothedDetector` — an exponential-moving-average wrapper that kills
  per-frame jitter; pure Python, so it (and the EMA maths) are testable without
  OpenCV.

Only the concrete detectors import ``cv2``; the smoothing layer does not.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

import numpy as np

from .models import FaceBox

YUNET_MODEL_URL = (
    "https://github.com/opencv/opencv_zoo/raw/main/models/"
    "face_detection_yunet/face_detection_yunet_2023mar.onnx"
)
YUNET_MODEL_NAME = "face_detection_yunet_2023mar.onnx"


class FaceDetector(Protocol):
    def detect(self, frame_bgr: np.ndarray) -> FaceBox | None: ...


def ema_box(prev: np.ndarray | None, target: np.ndarray, alpha: float) -> np.ndarray:
    """Exponential-moving-average blend of two [x, y, w, h] boxes."""
    if prev is None:
        return target
    return alpha * target + (1.0 - alpha) * prev


class SmoothedDetector:
    """Wrap any detector, smoothing its box over time to remove jitter."""

    def __init__(self, inner: FaceDetector, smoothing: float) -> None:
        self._inner = inner
        self._alpha = float(smoothing)
        self._state: np.ndarray | None = None

    def detect(self, frame_bgr: np.ndarray) -> FaceBox | None:
        box = self._inner.detect(frame_bgr)
        if box is None:
            return None
        target = np.array([box.x, box.y, box.w, box.h], dtype=np.float64)
        self._state = ema_box(self._state, target, self._alpha)
        x, y, w, h = self._state.astype(int)
        return FaceBox(int(x), int(y), int(w), int(h))


class HaarDetector:
    """Classic Haar-cascade detector (fallback)."""

    def __init__(self, cascade_path: Path) -> None:
        import cv2

        if not hasattr(cv2, "CascadeClassifier"):
            raise RuntimeError("this OpenCV build has no CascadeClassifier")
        self._cv2 = cv2
        self._cascade = cv2.CascadeClassifier(str(cascade_path))
        if self._cascade.empty():
            raise FileNotFoundError(f"Could not load Haar cascade: {cascade_path}")

    def detect(self, frame_bgr: np.ndarray) -> FaceBox | None:
        cv2 = self._cv2
        gray = cv2.cvtColor(frame_bgr, cv2.COLOR_BGR2GRAY)
        faces = self._cascade.detectMultiScale(
            gray, scaleFactor=1.1, minNeighbors=5, minSize=(60, 60)
        )
        if len(faces) == 0:
            return None
        x, y, w, h = max(faces, key=lambda f: f[2] * f[3])
        return FaceBox(int(x), int(y), int(w), int(h))


class YuNetDetector:
    """OpenCV YuNet DNN detector (preferred)."""

    def __init__(self, model_path: Path, score_threshold: float = 0.85) -> None:
        import cv2

        self._cv2 = cv2
        self._detector = cv2.FaceDetectorYN.create(
            str(model_path), "", (320, 320), score_threshold, 0.3, 5000
        )
        self._size: tuple[int, int] | None = None

    def detect(self, frame_bgr: np.ndarray) -> FaceBox | None:
        h, w = frame_bgr.shape[:2]
        if self._size != (w, h):
            self._detector.setInputSize((w, h))
            self._size = (w, h)
        _, faces = self._detector.detect(frame_bgr)
        if faces is None or len(faces) == 0:
            return None
        # Columns 0..3 are x, y, w, h; the largest is the person at the mirror.
        best = max(faces, key=lambda f: f[2] * f[3])
        x, y, bw, bh = best[:4]
        return FaceBox(int(x), int(y), int(bw), int(bh))


def ensure_yunet_model(cache_dir: Path, *, timeout: float = 30.0) -> Path:
    """Return a path to the YuNet model, downloading it once if needed."""
    import httpx

    model_dir = cache_dir / "models"
    model_dir.mkdir(parents=True, exist_ok=True)
    path = model_dir / YUNET_MODEL_NAME
    if path.exists() and path.stat().st_size > 0:
        return path
    response = httpx.get(YUNET_MODEL_URL, timeout=timeout, follow_redirects=True)
    response.raise_for_status()
    path.write_bytes(response.content)
    return path


def create_detector(settings, *, log=print) -> FaceDetector:
    """Build the best available detector, smoothed, per ``settings``.

    Tries YuNet (downloading its model on first use) and falls back to the Haar
    cascade if OpenCV lacks YuNet or the model can't be fetched.
    """
    inner: FaceDetector
    try:
        model = ensure_yunet_model(settings.cache_dir, timeout=settings.http_timeout)
        inner = YuNetDetector(model)
        log("Face detector: YuNet (DNN)")
    except Exception as yunet_exc:  # noqa: BLE001 - any failure => safe fallback
        try:
            inner = HaarDetector(settings.face_cascade)
            log(f"Face detector: Haar cascade (YuNet unavailable: {yunet_exc})")
        except Exception as haar_exc:  # noqa: BLE001
            raise RuntimeError(
                f"No face detector available (YuNet: {yunet_exc}; Haar: {haar_exc}). "
                "Install a full OpenCV build: pip install 'genetic-reflexions[camera]'."
            ) from haar_exc
    return SmoothedDetector(inner, settings.face_smoothing)
