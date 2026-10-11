from __future__ import annotations

import numpy as np

from magic_mirror.config import Settings
from magic_mirror.models import FaceBox, Trait
from magic_mirror.render import Renderer


def _traits():
    return [
        Trait("rs53576", "GG", 2.5, "Optimistic and empathetic", "cognition"),
        Trait("rs1815739", "CC", 2.2, "Better performing muscles", "performance"),
    ]


def test_renderer_render_preserves_shape_and_draws():
    renderer = Renderer(Settings(), header="Test", subheader="sub")
    frame = np.zeros((240, 320, 3), dtype=np.uint8)
    out = renderer.render(frame, _traits(), time_s=1.2, face=FaceBox(120, 80, 80, 90))
    assert out.shape == frame.shape
    assert out.dtype == np.uint8
    assert out.sum() > 0


def test_renderer_handles_no_traits_and_no_face():
    renderer = Renderer(Settings(), header="Test", subheader="sub")
    frame = np.zeros((200, 200, 3), dtype=np.uint8)
    out = renderer.render(frame, [], time_s=0.0, face=None)
    assert out.shape == frame.shape
