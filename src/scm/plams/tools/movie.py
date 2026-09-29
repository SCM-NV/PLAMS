import base64
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import (
    Iterable,
    Iterator,
    Optional,
    Sequence,
    Union,
    TYPE_CHECKING,
)

from scm.plams.core.functions import requires_optional_package

if TYPE_CHECKING:
    from PIL import Image as PilImage

__all__ = ["movie", "display_movie"]

Frames = Union[str, os.PathLike, Iterable[object]]

_IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".bmp", ".tif", ".tiff", ".webp"}


@requires_optional_package("PIL")
def movie(
    frames: Frames,
    output: Union[str, os.PathLike],
    *,
    fps: float = 10,
    loop: int = 0,
    sort: bool = True,
    ffmpeg: Optional[Union[str, os.PathLike]] = None,
    keep_frames: bool = False,
    matplotlib_dpi: int = 150,
) -> Path:
    """
    Write image frames to a GIF, animated WebP, or MP4 movie.

    ``frames`` can be a directory of images, an image file path, a sequence or generator
    of image file paths, PIL images, matplotlib figures, or matplotlib axes.

    :param frames: input frames
    :param output: output movie path, ending in ``.gif``, ``.webp``, or ``.mp4``
    :param fps: frames per second
    :param loop: GIF loop count, where ``0`` means loop forever
    :param sort: sort directory frames lexically
    :param ffmpeg: optional path to the ffmpeg executable for MP4 output
    :param keep_frames: keep PNG frames written next to the output movie
    :param matplotlib_dpi: DPI used when rendering matplotlib objects
    :return: path to the written movie
    """
    output_path = Path(output)
    suffix = output_path.suffix.lower()

    if fps <= 0:
        raise ValueError(f"fps must be positive, but was {fps!r}")
    if suffix not in {".gif", ".webp", ".mp4"}:
        raise ValueError(f"Unsupported movie format '{output_path.suffix}'. Supported formats are: .gif, .webp, .mp4")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if suffix in {".gif", ".webp"}:
        frame_dir = (
            Path(tempfile.mkdtemp(prefix=f"{output_path.stem}_frames_", dir=output_path.parent))
            if keep_frames
            else None
        )
        _write_pillow_movie(
            frames,
            output_path,
            fps=fps,
            loop=loop,
            sort=sort,
            frame_dir=frame_dir,
            matplotlib_dpi=matplotlib_dpi,
        )
    elif keep_frames:
        frame_dir = Path(tempfile.mkdtemp(prefix=f"{output_path.stem}_frames_", dir=output_path.parent))
        _write_mp4_from_frame_dir(
            _find_ffmpeg(ffmpeg),
            _write_png_frames(frames, frame_dir, sort=sort, matplotlib_dpi=matplotlib_dpi),
            output_path,
            fps=fps,
        )
    else:
        with tempfile.TemporaryDirectory(prefix="plams_movie_") as tmpdir:
            _write_mp4_from_frame_dir(
                _find_ffmpeg(ffmpeg),
                _write_png_frames(frames, Path(tmpdir), sort=sort, matplotlib_dpi=matplotlib_dpi),
                output_path,
                fps=fps,
            )

    return output_path


@requires_optional_package("IPython")
def display_movie(
    path: Union[str, os.PathLike],
    *,
    embed: bool = True,
    width: Optional[int] = None,
    height: Optional[int] = None,
) -> None:
    """
    Display a GIF, animated WebP, or MP4 movie in a Jupyter notebook.

    :param path: path to a GIF, WebP, or MP4 movie
    :param embed: embed MP4 video data in the notebook, defaults to ``True``
    :param width: optional display width in pixels
    :param height: optional display height in pixels
    """
    from IPython.display import HTML, Image, Video, display

    movie_path = Path(path)
    suffix = movie_path.suffix.lower()

    if suffix == ".gif":
        display(Image(filename=str(movie_path), width=width, height=height))
    elif suffix == ".webp":
        try:
            display(Image(filename=str(movie_path), width=width, height=height))
        except ValueError as exc:
            if "webp" not in str(exc).lower():
                raise
            data = base64.b64encode(movie_path.read_bytes()).decode("ascii")
            attributes = [f'src="data:image/webp;base64,{data}"']
            if width is not None:
                attributes.append(f'width="{width}"')
            if height is not None:
                attributes.append(f'height="{height}"')
            display(HTML(f"<img {' '.join(attributes)}>"))
    elif suffix == ".mp4":
        display(Video(str(movie_path), embed=embed, width=width, height=height))
    else:
        raise ValueError(f"Unsupported movie format '{movie_path.suffix}'. Supported formats are: .gif, .webp, .mp4")


def _iter_frames(frames: Frames, *, sort: bool) -> Iterator[object]:
    if isinstance(frames, (str, os.PathLike)):
        path = Path(frames)
        if path.is_dir():
            paths = [p for p in path.iterdir() if p.is_file() and p.suffix.lower() in _IMAGE_SUFFIXES]
            if sort:
                paths = sorted(paths)
            if not paths:
                raise ValueError(f"No supported image files found in directory: {path}")
            yield from paths
        else:
            yield path
        return

    try:
        yield from frames
    except TypeError as exc:
        raise TypeError(
            "frames must be a path, directory, or iterable of image paths, PIL images, or matplotlib objects"
        ) from exc


def _write_pillow_movie(
    frames: Frames,
    output: Path,
    *,
    fps: float,
    loop: int,
    sort: bool,
    frame_dir: Optional[Path],
    matplotlib_dpi: int,
) -> None:
    images: Sequence["PilImage.Image"]
    if frame_dir is not None:
        frame_paths = _write_png_frames(frames, frame_dir, sort=sort, matplotlib_dpi=matplotlib_dpi)
        images = [_frame_to_image(path, matplotlib_dpi=matplotlib_dpi) for path in frame_paths]
    else:
        images = [_frame_to_image(frame, matplotlib_dpi=matplotlib_dpi) for frame in _iter_frames(frames, sort=sort)]

    if not images:
        raise ValueError("No frames were supplied")

    images = _normalize_images_to_canvas(images)
    first, rest = images[0], images[1:]
    if output.suffix.lower() == ".webp":
        first.save(
            output,
            format="WEBP",
            save_all=True,
            append_images=rest,
            duration=round(1000 / fps),
            loop=loop,
        )
        return

    first.save(
        output,
        save_all=True,
        append_images=rest,
        duration=round(1000 / fps),
        loop=loop,
        disposal=2,
        background=0,
        optimize=False,
    )


def _write_mp4_from_frame_dir(
    ffmpeg_path: Path,
    frame_paths: Sequence[Path],
    output: Path,
    *,
    fps: float,
) -> None:
    if not frame_paths:
        raise ValueError("No frames were supplied")

    command = [
        str(ffmpeg_path),
        "-y",
        "-v",
        "error",
        "-nostdin",
        "-hide_banner",
        "-framerate",
        str(fps),
        "-i",
        str(frame_paths[0].parent / "frame_%06d.png"),
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-pix_fmt",
        "yuv420p",
        str(output),
    ]
    result = subprocess.run(
        command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        error = result.stderr.strip() if result.stderr else "unknown error"
        raise RuntimeError(f"ffmpeg failed while writing '{output}': {error}")


def _write_png_frames(frames: Frames, frame_dir: Path, *, sort: bool, matplotlib_dpi: int) -> Sequence[Path]:
    frame_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = []
    images = [_frame_to_image(frame, matplotlib_dpi=matplotlib_dpi) for frame in _iter_frames(frames, sort=sort)]
    for index, img in enumerate(_normalize_images_to_canvas(images)):
        path = frame_dir / f"frame_{index:06d}.png"
        img.save(path)
        frame_paths.append(path)
    return frame_paths


def _frame_to_image(frame: object, *, matplotlib_dpi: int) -> "PilImage.Image":
    from PIL import Image as PilImage

    if isinstance(frame, PilImage.Image):
        return frame.copy()

    if isinstance(frame, (str, os.PathLike)):
        path = Path(frame)
        if not path.is_file():
            raise FileNotFoundError(f"Frame image does not exist: {path}")
        with PilImage.open(path) as img:
            return img.copy()

    figure = frame if hasattr(frame, "savefig") and hasattr(frame, "canvas") else getattr(frame, "figure", None)
    if figure is not None and hasattr(figure, "savefig") and hasattr(figure, "canvas"):
        return _matplotlib_frame_to_image(frame, dpi=matplotlib_dpi)

    raise TypeError(
        "Frame must be an image path, PIL image, matplotlib figure, or matplotlib axes; " f"got {type(frame).__name__}"
    )


def _as_rgb_image(img: "PilImage.Image") -> "PilImage.Image":
    from PIL import Image as PilImage

    if img.mode == "RGB":
        return img

    if img.mode == "RGBA" or (img.mode == "P" and "transparency" in img.info):
        rgba = img.convert("RGBA")
        background = PilImage.new("RGBA", rgba.size, "white")
        background.alpha_composite(rgba)
        return background.convert("RGB")

    return img.convert("RGB")


def _normalize_images_to_canvas(
    images: Sequence["PilImage.Image"],
) -> Sequence["PilImage.Image"]:
    from PIL import Image as PilImage

    rgb_images = [_as_rgb_image(img) for img in images]
    if not rgb_images:
        return rgb_images

    width = max(img.width for img in rgb_images)
    height = max(img.height for img in rgb_images)
    normalized = []
    for img in rgb_images:
        if img.size == (width, height):
            normalized.append(img)
            continue

        canvas = PilImage.new("RGB", (width, height), "white")
        canvas.paste(img, ((width - img.width) // 2, (height - img.height) // 2))
        normalized.append(canvas)
    return normalized


@requires_optional_package("matplotlib")
def _matplotlib_frame_to_image(frame: object, *, dpi: int) -> "PilImage.Image":
    from io import BytesIO

    from PIL import Image as PilImage

    figure = frame if hasattr(frame, "savefig") and hasattr(frame, "canvas") else getattr(frame, "figure")
    buffer = BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi)
    buffer.seek(0)
    with PilImage.open(buffer) as img:
        return img.copy()


def _find_ffmpeg(ffmpeg: Optional[Union[str, os.PathLike]]) -> Path:
    candidates = []
    if ffmpeg is not None:
        candidates.append(Path(ffmpeg))
    if os.environ.get("SCM_FFMPEG"):
        candidates.append(Path(os.environ["SCM_FFMPEG"]))
    if os.environ.get("AMSBIN"):
        candidates.extend(
            [
                Path(os.environ["AMSBIN"]) / "ffmpeg.exe",
                Path(os.environ["AMSBIN"]) / "ffmpeg",
            ]
        )
    candidates.extend(
        [
            Path("ffmpeg"),
            Path("ffmpeg.exe"),
            Path("/usr/local/bin/ffmpeg"),
            Path("/usr/local/bin/ffmpeg.exe"),
        ]
    )

    for candidate in candidates:
        if candidate.is_file():
            return candidate
        if not candidate.is_absolute() and candidate.parent == Path("."):
            path = shutil.which(os.fspath(candidate))
            if path:
                return Path(path)

    raise RuntimeError("Writing MP4 movies requires ffmpeg. Set SCM_FFMPEG, pass ffmpeg=..., or write a GIF instead.")
