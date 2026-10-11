from __future__ import annotations

import numpy as np
import pytest

from magic_mirror import detect
from magic_mirror.config import Settings
from magic_mirror.detect import SmoothedDetector, ema_box
from magic_mirror.models import FaceBox


def test_ema_box_first_reading_passes_through():
    target = np.array([10, 20, 30, 40], dtype=float)
    assert np.array_equal(ema_box(None, target, 0.5), target)


def test_ema_box_blends():
    prev = np.array([0, 0, 0, 0], dtype=float)
    target = np.array([10, 10, 10, 10], dtype=float)
    assert np.allclose(ema_box(prev, target, 0.5), [5, 5, 5, 5])


class _FixedDetector:
    def __init__(self, box):
        self.box = box

    def detect(self, _frame):
        return self.box


def test_smoothed_detector_converges():
    inner = _FixedDetector(FaceBox(100, 100, 50, 50))
    smoothed = SmoothedDetector(inner, smoothing=0.5)
    frame = np.zeros((10, 10, 3), dtype=np.uint8)

    first = smoothed.detect(frame)
    assert first == FaceBox(100, 100, 50, 50)  # seeds on first reading
    # Repeated identical readings stay put (already converged).
    assert smoothed.detect(frame) == FaceBox(100, 100, 50, 50)


def test_smoothed_detector_passes_through_none():
    smoothed = SmoothedDetector(_FixedDetector(None), smoothing=0.5)
    assert smoothed.detect(np.zeros((10, 10, 3), dtype=np.uint8)) is None


def test_create_detector_falls_back_to_haar(tmp_path, monkeypatch):
    cv2 = pytest.importorskip("cv2")  # Haar detector needs OpenCV
    if not hasattr(cv2, "CascadeClassifier"):
        pytest.skip("this OpenCV build has no Haar CascadeClassifier")

    def unavailable(*_a, **_k):
        raise RuntimeError("no model")

    monkeypatch.setattr(detect, "ensure_yunet_model", unavailable)
    settings = Settings(cache_dir=tmp_path)  # face_cascade defaults to a real file
    detector = detect.create_detector(settings, log=lambda *_: None)

    assert isinstance(detector, SmoothedDetector)
    assert isinstance(detector._inner, detect.HaarDetector)
    # Runs on a blank frame without raising (finds nothing).
    assert detector.detect(np.zeros((80, 80, 3), dtype=np.uint8)) is None
