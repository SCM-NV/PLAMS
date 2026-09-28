#!/usr/bin/env amspython
# coding: utf-8

import os
import subprocess
import base64
from pathlib import Path
from unittest.mock import patch

import pytest
from PIL import Image as PilImage

from scm.plams.tools.movie import display_movie, movie


def make_image(color):
    return PilImage.new("RGB", (8, 6), color=color)


def test_movie_functions_are_exported_from_plams():
    import scm.plams as plams

    assert plams.movie is movie
    assert plams.display_movie is display_movie


def test_movie_writes_gif_from_pil_images(tmp_path):
    output = tmp_path / "movie.gif"

    result = movie([make_image("red"), make_image("blue")], output, fps=5)

    assert result == output
    with PilImage.open(output) as img:
        assert img.format == "GIF"
        assert img.n_frames == 2


def test_movie_writes_animated_webp_from_pil_images(tmp_path):
    output = tmp_path / "movie.webp"

    result = movie([make_image("red"), make_image("blue")], output, fps=5)

    assert result == output
    with PilImage.open(output) as img:
        assert img.format == "WEBP"
        assert img.n_frames == 2


def test_movie_writes_gif_from_directory_in_sorted_order(tmp_path):
    make_image("blue").save(tmp_path / "frame_002.png")
    make_image("red").save(tmp_path / "frame_001.png")
    output = tmp_path / "movie.gif"

    movie(tmp_path, output)

    with PilImage.open(output) as img:
        img.seek(0)
        assert img.convert("RGB").getpixel((0, 0)) == (255, 0, 0)
        img.seek(1)
        assert img.convert("RGB").getpixel((0, 0)) == (0, 0, 255)


def test_movie_writes_full_gif_frames_for_transparent_images(tmp_path):
    red = PilImage.new("RGBA", (8, 6), color=(255, 255, 255, 0))
    red.putpixel((0, 0), (255, 0, 0, 255))
    blue = PilImage.new("RGBA", (8, 6), color=(255, 255, 255, 0))
    blue.putpixel((7, 0), (0, 0, 255, 255))
    output = tmp_path / "movie.gif"

    movie([red, blue], output)

    with PilImage.open(output) as img:
        img.seek(1)
        frame = img.convert("RGB")
        assert frame.getpixel((0, 0)) == (255, 255, 255)
        assert frame.getpixel((7, 0)) == (0, 0, 255)


def test_movie_writes_gif_frames_on_same_full_canvas(tmp_path):
    small = PilImage.new("RGBA", (6, 4), color=(255, 255, 255, 0))
    small.putpixel((0, 0), (255, 0, 0, 255))
    large = PilImage.new("RGBA", (10, 8), color=(255, 255, 255, 0))
    large.putpixel((9, 7), (0, 0, 255, 255))
    output = tmp_path / "movie.gif"

    movie([small, large], output)

    with PilImage.open(output) as img:
        assert img.size == (10, 8)
        for i in range(img.n_frames):
            img.seek(i)
            assert img.tile[0][1] == (0, 0, 10, 8)
        img.seek(0)
        frame = img.convert("RGB")
        assert frame.getpixel((2, 2)) == (255, 0, 0)
        assert frame.getpixel((0, 0)) == (255, 255, 255)


def test_movie_accepts_generator_of_image_paths(tmp_path):
    frame_paths = []
    for index, color in enumerate(["red", "green", "blue"]):
        path = tmp_path / f"frame_{index}.png"
        make_image(color).save(path)
        frame_paths.append(path)
    output = tmp_path / "movie.gif"

    movie((path for path in frame_paths), output)

    with PilImage.open(output) as img:
        assert img.n_frames == 3


def test_movie_accepts_matplotlib_axes(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    output = tmp_path / "movie.gif"

    try:
        movie([ax], output)
    finally:
        plt.close(fig)

    with PilImage.open(output) as img:
        assert img.format == "GIF"


def test_movie_rejects_empty_frames(tmp_path):
    with pytest.raises(ValueError, match="No frames"):
        movie([], tmp_path / "movie.gif")


def test_movie_rejects_unsupported_output_suffix(tmp_path):
    with pytest.raises(ValueError, match="Unsupported movie format"):
        movie([make_image("red")], tmp_path / "movie.avi")


def test_movie_mp4_requires_ffmpeg(tmp_path, monkeypatch):
    monkeypatch.delenv("SCM_FFMPEG", raising=False)
    monkeypatch.delenv("AMSBIN", raising=False)
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(RuntimeError, match="requires ffmpeg"):
        movie([make_image("red")], tmp_path / "movie.mp4")


def test_movie_mp4_runs_ffmpeg_with_png_frames(tmp_path, monkeypatch):
    fake_ffmpeg = tmp_path / "ffmpeg"
    fake_ffmpeg.write_text("", encoding="utf-8")
    output = tmp_path / "movie.mp4"

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        movie([make_image("red"), make_image("blue")], output, fps=12, ffmpeg=fake_ffmpeg)

    command = mock_run.call_args.args[0]
    assert command[0] == os.fspath(fake_ffmpeg)
    assert command[1:6] == [
        "-y",
        "-v",
        "error",
        "-nostdin",
        "-hide_banner",
    ]
    assert "-framerate" in command
    assert command[command.index("-framerate") + 1] == "12"
    assert command[-1] == os.fspath(output)
    mock_run.assert_called_once()
    assert mock_run.call_args.kwargs == {
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.PIPE,
        "text": True,
    }


def test_movie_mp4_can_keep_frame_directory(tmp_path):
    fake_ffmpeg = tmp_path / "ffmpeg"
    fake_ffmpeg.write_text("", encoding="utf-8")

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 0
        movie(
            [make_image("red")],
            tmp_path / "movie.mp4",
            ffmpeg=fake_ffmpeg,
            keep_frames=True,
        )

    frame_dirs = list(Path(tmp_path).glob("movie_frames_*"))
    assert len(frame_dirs) == 1
    assert (frame_dirs[0] / "frame_000000.png").is_file()


def test_movie_mp4_reports_ffmpeg_errors(tmp_path):
    fake_ffmpeg = tmp_path / "ffmpeg"
    fake_ffmpeg.write_text("", encoding="utf-8")

    with patch("subprocess.run") as mock_run:
        mock_run.return_value.returncode = 1
        mock_run.return_value.stderr = "bad codec"
        with pytest.raises(RuntimeError, match="bad codec"):
            movie([make_image("red")], tmp_path / "movie.mp4", ffmpeg=fake_ffmpeg)


def test_movie_gif_can_keep_frame_directory(tmp_path):
    movie([make_image("red")], tmp_path / "movie.gif", keep_frames=True)

    frame_dirs = list(Path(tmp_path).glob("movie_frames_*"))
    assert len(frame_dirs) == 1
    assert (frame_dirs[0] / "frame_000000.png").is_file()


def test_display_movie_displays_gif(tmp_path):
    pytest.importorskip("IPython")
    path = tmp_path / "movie.gif"
    movie([make_image("red")], path)

    with patch("IPython.display.display") as mock_display:
        result = display_movie(path, width=320, height=240)

    assert result is None
    displayed = mock_display.call_args.args[0]
    assert displayed.filename == os.fspath(path)
    assert displayed.width == 320
    assert displayed.height == 240


def test_display_movie_displays_webp(tmp_path):
    pytest.importorskip("IPython")
    path = tmp_path / "movie.webp"
    path.write_bytes(b"fake webp")

    with patch("IPython.display.Image") as mock_image, patch("IPython.display.display") as mock_display:
        result = display_movie(path, width=320)

    assert result is None
    mock_image.assert_called_once_with(filename=os.fspath(path), width=320, height=None)
    assert mock_display.call_args.args[0] is mock_image.return_value


def test_display_movie_falls_back_to_html_for_old_ipython_webp(tmp_path):
    pytest.importorskip("IPython")
    path = tmp_path / "movie.webp"
    data = b"fake webp"
    path.write_bytes(data)

    with patch(
        "IPython.display.Image",
        side_effect=ValueError("Cannot embed the 'webp' image format"),
    ), patch("IPython.display.display") as mock_display:
        result = display_movie(path, width=320)

    assert result is None
    displayed = mock_display.call_args.args[0]
    expected = base64.b64encode(data).decode("ascii")
    assert f'src="data:image/webp;base64,{expected}"' in displayed.data
    assert 'width="320"' in displayed.data


def test_display_movie_displays_mp4(tmp_path):
    pytest.importorskip("IPython")
    path = tmp_path / "movie.mp4"
    path.write_bytes(b"fake mp4")

    with patch("IPython.display.display") as mock_display:
        result = display_movie(path, embed=False, width=640)

    assert result is None
    displayed = mock_display.call_args.args[0]
    assert displayed.filename == os.fspath(path)
    assert displayed.embed is False
    assert displayed.width == 640


def test_display_movie_rejects_unsupported_suffix(tmp_path):
    pytest.importorskip("IPython")

    with pytest.raises(ValueError, match="Unsupported movie format"):
        display_movie(tmp_path / "movie.avi")
