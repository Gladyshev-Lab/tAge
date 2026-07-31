"""Generate the tAge logo, as SVG only.

Two assets come out of this, because they are consumed at very different sizes:

    man/figures/logo-full.svg   clock + wordmark. Legible at README size, which
                                is its one and only use.
    man/figures/logo.svg        square clock mark. pkgdown's find_logo() prefers
                                logo.svg over logo.png, so this is what the
                                navbar shows, and it is the source the favicon
                                set is generated from. The wordmark is not in it
                                because it is unreadable below ~48 px.

Nothing here rasterises: the geometry is emitted straight to SVG, and the mark's
viewBox is computed from the drawn extents rather than measured off a bitmap.
Pillow is still imported, but only for font metrics - see SvgCanvas.text().

After changing the artwork, regenerate the icon set (7 files, never hand-edited):

    make favicons
"""

import os
import shutil

import numpy as np
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.ttLib import TTFont
from PIL import ImageFont

# -- Canvas ------------------------------------------------------------------

width, height = 600, 600

# Color scheme
primary_color = (30, 100, 180)
accent_color = (100, 180, 230)

center_x, center_y = width // 2, height // 2
center_y -= 50

num_points = 50
clock_radius = 100
clock_y = center_y - 50


class SvgCanvas:
    """Collects drawing primitives as SVG elements, tracking their extents.

    The method signatures mirror the subset of PIL's ImageDraw that the artwork
    uses, so draw_helix() and draw_dial() describe the geometry exactly once and
    know nothing about the output format.
    """

    def __init__(self):
        self.parts = []
        self.bounds = None      # (min_x, min_y, max_x, max_y)

    def _grow(self, x0, y0, x1, y1):
        if self.bounds is None:
            self.bounds = [x0, y0, x1, y1]
        else:
            b = self.bounds
            b[0], b[1] = min(b[0], x0), min(b[1], y0)
            b[2], b[3] = max(b[2], x1), max(b[3], y1)

    @staticmethod
    def _paint(color):
        r, g, b = color[:3]
        alpha = color[3] / 255 if len(color) > 3 else 1.0
        opacity = '' if alpha >= 1 else f' opacity="{alpha:.3f}"'
        return f'rgb({r},{g},{b})', opacity

    def ellipse(self, box, fill=None, outline=None, width=1):
        x0, y0, x1, y1 = box
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        rx, ry = (x1 - x0) / 2, (y1 - y0) / 2
        if outline is not None:
            # An SVG stroke straddles the path, so pull the radius in by half
            # the stroke width to keep the outer edge on the given box.
            rx -= width / 2
            ry -= width / 2
        attrs = f'cx="{cx:.2f}" cy="{cy:.2f}" rx="{rx:.2f}" ry="{ry:.2f}"'
        if fill is not None:
            colour, opacity = self._paint(fill)
            attrs += f' fill="{colour}"{opacity}'
        else:
            attrs += ' fill="none"'
        if outline is not None:
            colour, opacity = self._paint(outline)
            attrs += f' stroke="{colour}" stroke-width="{width}"{opacity}'
        self.parts.append(f'<ellipse {attrs}/>')
        self._grow(x0, y0, x1, y1)

    def line(self, points, fill=None, width=1):
        x0, y0, x1, y1 = points
        colour, opacity = self._paint(fill)
        self.parts.append(
            f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" '
            f'stroke="{colour}" stroke-width="{width}" '
            f'stroke-linecap="round"{opacity}/>'
        )
        half = width / 2      # round caps stick out in every direction
        self._grow(min(x0, x1) - half, min(y0, y1) - half,
                   max(x0, x1) + half, max(y0, y1) + half)

    def text(self, xy, string, font, fill):
        """Add a string as outlines rather than an SVG <text> element.

        <text> would depend on DejaVu Sans being installed wherever the logo is
        viewed and silently fall back to something else where it is not. Glyph
        outlines render identically everywhere.

        Pen positions come from Pillow's getlength(), which applies the font's
        kerning; stepping through raw advance widths would not.
        """
        tt = TTFont(font.path, fontNumber=getattr(font, 'index', 0))
        upem = tt['head'].unitsPerEm
        scale = font.size / upem
        ascent, _ = font.getmetrics()
        baseline = xy[1] + ascent
        cmap = tt.getBestCmap()
        glyphs = tt.getGlyphSet()
        colour, opacity = self._paint(fill)

        for i, char in enumerate(string):
            name = cmap.get(ord(char))
            if name is None:
                continue
            pen = SVGPathPen(glyphs)
            glyphs[name].draw(pen)
            commands = pen.getCommands()
            if not commands:          # spaces and other blank glyphs
                continue
            x = xy[0] + font.getlength(string[:i])
            self.parts.append(
                f'<path d="{commands}" fill="{colour}"{opacity} '
                f'transform="translate({x:.2f} {baseline:.2f}) '
                f'scale({scale:.6f} {-scale:.6f})"/>'
            )
        tt.close()

    def save(self, path, viewbox):
        x, y, w, h = viewbox
        body = '\n  '.join(self.parts)
        with open(path, 'w', encoding='utf-8') as fh:
            fh.write(
                f'<svg xmlns="http://www.w3.org/2000/svg" '
                f'viewBox="{x:g} {y:g} {w:g} {h:g}" '
                f'width="{w:g}" height="{h:g}">\n'
                f'  {body}\n</svg>\n'
            )


# -- Artwork -----------------------------------------------------------------


def draw_helix(target, cx, cy, r_start, r_step, dot_base=5, dot_growth=0.06,
               points=num_points):
    """Draw the two-strand spiral of dots that rings the dial.

    Parameterised because the small icon needs the same motif compressed: at
    favicon sizes a ring sitting far out from the dial thins to nothing, so the
    mark redraws it close in rather than reusing the large version's radii.
    """
    for i in range(points):
        angle = (i / points) * 3 * np.pi
        radius = r_start + i * r_step
        x1 = cx + radius * np.cos(angle)
        y1 = cy + radius * np.sin(angle)
        x2 = cx + radius * np.cos(angle + np.pi)
        y2 = cy + radius * np.sin(angle + np.pi)

        size = int(dot_base + i * dot_growth)
        alpha = int(180 + (i / points) * 75)

        target.ellipse([x1-size, y1-size, x1+size, y1+size],
                       fill=primary_color + (alpha,))
        target.ellipse([x2-size, y2-size, x2+size, y2+size],
                       fill=accent_color + (alpha,))


def draw_dial(target, cx, cy, radius=clock_radius):
    """Draw the dial: circle, two hands, centre dot."""
    target.ellipse([cx - radius, cy - radius, cx + radius, cy + radius],
                   outline=primary_color + (255,), width=6)

    hand_angle = -np.pi / 3
    hand_length = 60
    target.line([cx, cy,
                 cx + hand_length * np.cos(hand_angle),
                 cy + hand_length * np.sin(hand_angle)],
                fill=primary_color + (255,), width=5)

    minute_angle = np.pi / 6
    minute_length = 80
    target.line([cx, cy,
                 cx + minute_length * np.cos(minute_angle),
                 cy + minute_length * np.sin(minute_angle)],
                fill=accent_color + (255,), width=4)

    target.ellipse([cx - 8, cy - 8, cx + 8, cy + 8],
                   fill=primary_color + (255,))


# Fonts. The Debian-style paths do not exist on Arch, but Pillow falls back to
# searching the system font directories by basename, so both resolve.
try:
    font_large = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 120)
    font_small = ImageFont.truetype(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 36)
except OSError as exc:
    raise SystemExit(f"DejaVu Sans not found, cannot draw the wordmark: {exc}")

text = "tAge"
text_x = (width - font_large.getlength(text)) // 2
text_y = height - 200

subtitle = "Transcriptomic Age"
sub_x = (width - font_small.getlength(subtitle)) // 2
sub_y = text_y + 130


# -- Output ------------------------------------------------------------------

FULL_PATH = 'man/figures/logo-full.svg'
MARK_PATH = 'man/figures/logo.svg'
# The Python docs use the same two files; kept in sync here so the sites cannot
# drift apart.
PYTHON_STATIC = '../tage-python/docs/_static'

MARK_PADDING = 12

# What the mark contains:
#   'compact'  dial plus the ring of dots, redrawn close in. Keeps the motif of
#              the full logo and still reads at 16 px, because the dial fills
#              most of the frame instead of floating inside a wide halo.
#   'wide'     dial plus the ring at the full logo's radii. Faithful to
#              logo-full.svg, but the outer dots thin to nothing and the whole
#              icon greys into a blob in a browser tab.
#   'clock'    dial only, no dots.
MARK_STYLE = 'compact'

# Radial band the dots occupy in the compact mark. The dial's outer edge is at
# 103 (radius 100 + half of the 6 px stroke), so the ring starts just clear of it
# and stops well short of where the full logo puts it (120 -> 179).
COMPACT_RING = dict(r_start=113, r_step=0.62, dot_base=5, dot_growth=0.05)
WIDE_RING = dict(r_start=120, r_step=1.2)

# Order matters: SVG has no z-index, it paints in document order.
full = SvgCanvas()
draw_helix(full, center_x, center_y - 50, **WIDE_RING)
draw_dial(full, center_x, clock_y)
full.text((text_x, text_y), text, font_large, primary_color + (255,))
full.text((sub_x, sub_y), subtitle, font_small, accent_color + (255,))
full.save(FULL_PATH, viewbox=(0, 0, width, height))
print(f"full logo  -> {FULL_PATH}  ({width}x{height})")

mark = SvgCanvas()
if MARK_STYLE == 'compact':
    draw_helix(mark, center_x, clock_y, **COMPACT_RING)
    draw_dial(mark, center_x, clock_y)
elif MARK_STYLE == 'clock':
    draw_dial(mark, center_x, clock_y)
elif MARK_STYLE == 'wide':
    draw_helix(mark, center_x, center_y - 50, **WIDE_RING)
    draw_dial(mark, center_x, clock_y)
else:
    raise SystemExit(f"unknown MARK_STYLE {MARK_STYLE!r}")

if mark.bounds is None:
    raise SystemExit(f"nothing drawn for MARK_STYLE={MARK_STYLE!r}")

# Square viewBox around whatever was drawn, so the mark is centred and tight
# regardless of which style produced it.
min_x, min_y, max_x, max_y = mark.bounds
mid_x, mid_y = (min_x + max_x) / 2, (min_y + max_y) / 2
side = max(max_x - min_x, max_y - min_y) + 2 * MARK_PADDING
mark.save(MARK_PATH, viewbox=(round(mid_x - side / 2, 2),
                              round(mid_y - side / 2, 2),
                              round(side, 2), round(side, 2)))
print(f"clock mark -> {MARK_PATH}  ({side:.0f}x{side:.0f}, style={MARK_STYLE})")

if os.path.isdir(PYTHON_STATIC):
    for src in (FULL_PATH, MARK_PATH):
        dst = os.path.join(PYTHON_STATIC, os.path.basename(src))
        shutil.copyfile(src, dst)
        print(f"           -> {dst}")
else:
    print(f"skipped {PYTHON_STATIC} (directory not found)")

print("\nNext, if the mark changed: make favicons")
print("  (the icon set stays raster - .ico and apple-touch-icon cannot be SVG)")
