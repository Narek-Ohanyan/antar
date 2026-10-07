"""Draws the contour-line artwork of the site from the project's own terrain raster (ui/assets/map/elevation.png, SRTM-derived, 650 m pixels).

Writes ui/assets/relief.svg: the national outline, minor contours every 250 m and major contours every 1000 m, in the map's own pixel space and cropped to
the country. Stroke colours are plain black with opacity, because the file is also used as a CSS mask (the page decoration takes its colour from the
theme); inline in the home page the same elements are recoloured by CSS. Nothing here is a model result: it is the terrain the models are run on.

    python3 scripts/build_relief_art.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402

UI = Path(__file__).resolve().parent.parent / "ui"
OUT = UI / "assets" / "relief.svg"
STEP_M = 250
MAJOR_M = 1000
SMOOTH = 2          # block-average the raster 2x2 before contouring: a calmer line at display size


def load():
    e = np.asarray(Image.open(UI / "assets" / "map" / "elevation.png").convert("RGB"), dtype=float)
    z = e[..., 0] * 256 + e[..., 1]
    inside = np.asarray(Image.open(UI / "assets" / "map" / "region.png").convert("L")) > 0
    return z, inside


def block_mean(a, k):
    h, w = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    return a[:h, :w].reshape(h // k, k, w // k, k).mean(axis=(1, 3))


def rdp(xy, eps):
    """Douglas-Peucker line simplification: keeps every point that lies more than ``eps`` pixels off the chord."""
    if len(xy) < 3:
        return xy
    a, b = xy[0], xy[-1]
    ab = b - a
    n = np.hypot(*ab)
    d = np.abs(ab[0] * (xy[:, 1] - a[1]) - ab[1] * (xy[:, 0] - a[0])) / n if n > 0 else np.hypot(*(xy - a).T)
    i = int(np.argmax(d))
    if d[i] <= eps:
        return np.array([a, b])
    left, right = rdp(xy[: i + 1], eps), rdp(xy[i:], eps)
    return np.vstack([left[:-1], right])


def paths_d(cs, eps=0.45):
    parts = []
    for seg_list in cs.allsegs:
        for seg in seg_list:
            if len(seg) < 4:
                continue
            xy = rdp(np.asarray(seg), eps)
            if len(xy) < 2:
                continue
            parts.append("M" + "L".join(f"{x:.1f} {y:.1f}" for x, y in xy))
    return "".join(parts)


def build():
    z, inside = load()
    h, w = z.shape
    bbox_y, bbox_x = np.where(inside)
    pad = 6
    x0, x1, y0, y1 = max(bbox_x.min() - pad, 0), min(bbox_x.max() + pad, w), max(bbox_y.min() - pad, 0), min(bbox_y.max() + pad, h)

    zin = np.where(inside, z, np.nan)
    frac = block_mean(inside.astype(float), SMOOTH)
    zs = block_mean(np.where(inside, z, 0.0), SMOOTH) / np.maximum(frac, 1e-9)
    zs = np.ma.masked_where(frac < 0.75, zs)                      # a smoothed pixel is kept only if at least 3 of its 4 source pixels are in the country
    gy, gx = np.mgrid[0:zs.shape[0], 0:zs.shape[1]]
    X, Y = (gx + 0.5) * SMOOTH, (gy + 0.5) * SMOOTH

    fig, ax = plt.subplots()
    levels = np.arange(STEP_M, np.nanmax(zin) + STEP_M, STEP_M)
    minor = [v for v in levels if v % MAJOR_M]
    major = [v for v in levels if v % MAJOR_M == 0]
    cs_minor = ax.contour(X, Y, zs, levels=minor)
    cs_major = ax.contour(X, Y, zs, levels=major)
    cs_out = ax.contour(np.arange(w) + 0.5, np.arange(h) + 0.5, inside.astype(float), levels=[0.5])
    d_minor, d_major, d_out = paths_d(cs_minor), paths_d(cs_major), paths_d(cs_out)
    plt.close(fig)

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{x0} {y0} {x1 - x0} {y1 - y0}" fill="none" stroke="#000" stroke-linecap="round" stroke-linejoin="round">'
           f'<title>Contour lines of Armenia every {STEP_M} m</title>'
           f'<path class="minor" d="{d_minor}" stroke-width="0.7" stroke-opacity="0.55"/>'
           f'<path class="major" d="{d_major}" stroke-width="1.1" stroke-opacity="0.9"/>'
           f'<path class="outline" d="{d_out}" stroke-width="1.6"/></svg>\n')
    OUT.write_text(svg)
    print(f"wrote {OUT} ({OUT.stat().st_size / 1024:.0f} kB); viewBox {x0} {y0} {x1 - x0} {y1 - y0}; levels {len(minor)} minor, {len(major)} major")


if __name__ == "__main__":
    build()
