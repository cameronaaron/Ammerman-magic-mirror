from __future__ import annotations

import re

from typer.testing import CliRunner

from magic_mirror.cli import app

runner = CliRunner()


def test_list_csv_source():
    result = runner.invoke(app, ["--source", "csv", "--no-usb", "list"])
    assert result.exit_code == 0, result.output
    # The bundled sample is synthetic, so assert that notable markers are
    # listed, not which ones (the old check named a marker of a real genome).
    assert re.search(r"\brs\d+\b", result.output)
    assert "notable" in result.output


def test_poster_csv_source(tmp_path):
    out = tmp_path / "poster.png"
    result = runner.invoke(
        app,
        ["--source", "csv", "--no-usb", "poster", str(out), "--width", "640", "--height", "400"],
    )
    assert result.exit_code == 0, result.output
    assert out.exists() and out.stat().st_size > 0


def test_list_bad_genome_errors(tmp_path):
    result = runner.invoke(
        app, ["--genome", str(tmp_path / "nope.txt"), "--no-usb", "--source", "csv", "list"]
    )
    assert result.exit_code == 1
    assert "error" in result.output.lower()


def test_devices_runs():
    result = runner.invoke(app, ["devices"])
    assert result.exit_code == 0


def test_help():
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "magic mirror" in result.output.lower()
