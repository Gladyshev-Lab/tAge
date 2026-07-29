from PIL import Image, ImageDraw, ImageFont
import numpy as np

# Create image with transparent background
width, height = 600, 600
img = Image.new('RGBA', (width, height), color=(255, 255, 255, 0))
draw = ImageDraw.Draw(img)

# Color scheme
primary_color = (30, 100, 180)
accent_color = (100, 180, 230)

center_x, center_y = width // 2, height // 2
center_y -= 50

# Simplified DNA helix - fewer points, cleaner look
num_points = 50


def draw_helix(target, cx, cy, r_start, r_step, dot_base=5, dot_growth=0.06,
               points=num_points):
    """Draw the two-strand spiral of dots that rings the dial.

    Parameterised because the small icon needs the same motif compressed: at
    favicon sizes a ring sitting far out from the dial thins to nothing, so the
    mark redraws it close in rather than cropping the large version.
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


# The helix is centred 50 px above the canvas centre, same as the dial below.
draw_helix(draw, center_x, center_y - 50, r_start=120, r_step=1.2)

# Simplified clock - cleaner design
clock_radius = 100
clock_y = center_y - 50


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


# The dial goes on its own transparent layer so the 'clock' mark can be taken
# from it cleanly; cropping the flattened image would drag fragments of the
# surrounding helix dots into the corners.
clock_layer = Image.new('RGBA', (width, height), color=(255, 255, 255, 0))
draw_dial(ImageDraw.Draw(clock_layer), center_x, clock_y)

# Flatten the dial onto the helix and go back to drawing on the main canvas.
img.alpha_composite(clock_layer)
draw = ImageDraw.Draw(img)

# Text
try:
    font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 120)
    font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 36)
except:
    font_large = ImageFont.load_default()
    font_small = ImageFont.load_default()

# Main text
text = "tAge"
bbox = draw.textbbox((0, 0), text, font=font_large)
text_width = bbox[2] - bbox[0]
text_x = (width - text_width) // 2
text_y = height - 200

draw.text((text_x, text_y), text, fill=primary_color + (255,), font=font_large)

# Subtitle
subtitle = "Transcriptomic Age"
bbox_sub = draw.textbbox((0, 0), subtitle, font=font_small)
sub_width = bbox_sub[2] - bbox_sub[0]
sub_x = (width - sub_width) // 2
sub_y = text_y + 130

draw.text((sub_x, sub_y), subtitle, fill=accent_color + (255,), font=font_small)

# ---------------------------------------------------------------------------
# Output: two assets, because they are consumed at very different sizes.
#
#   man/figures/logo-full.png   clock + wordmark. Legible only at README size
#                               (300 px), which is its one and only use.
#   man/figures/logo.png        square crop of the clock alone. pkgdown derives
#                               BOTH the navbar logo (~30 px) and every favicon
#                               from this file, and the wordmark is unreadable
#                               mush at those sizes.
#
# Neither crop is hardcoded: both are measured from the alpha channel, so if the
# artwork moves or changes size the crops follow. MARK_STYLE below picks what the
# small icon contains.
# ---------------------------------------------------------------------------
import os
import shutil

FULL_PATH = 'man/figures/logo-full.png'
MARK_PATH = 'man/figures/logo.png'
# The Python docs use the same mark for html_logo and html_favicon. Kept in sync
# here so the two sites cannot drift apart.
PYTHON_MARK_PATH = '../tage-python/docs/_static/logo.png'
MARK_PADDING = 12

# What the mark contains:
#   'compact'  dial plus the ring of dots, redrawn close in. Keeps the motif of
#              the full logo and still reads at 16 px, because the dial fills
#              most of the frame instead of floating inside a wide halo.
#   'wide'     the full artwork above the wordmark, cropped as-is. Faithful to
#              logo-full.png, but the outer dots thin to nothing and the whole
#              icon greys into a blob in a browser tab.
#   'clock'    dial only, no dots.
MARK_STYLE = 'compact'

# Radial band the dots occupy in the compact mark. The dial's outer edge is at
# 103 (radius 100 + half of the 6 px stroke), so the ring starts just clear of it
# and stops well short of where the full logo puts it (120 -> 179).
COMPACT_RING = dict(r_start=113, r_step=0.62, dot_base=5, dot_growth=0.05)

img.save(FULL_PATH)
print(f"full logo  -> {FULL_PATH}  ({img.width}x{img.height})")

if MARK_STYLE == 'compact':
    # Drawn fresh on its own square canvas rather than cropped: the ring has to
    # move inwards, which no crop of the large artwork can do.
    side_guess = 2 * (clock_radius + 60)
    source = Image.new('RGBA', (side_guess, side_guess), color=(255, 255, 255, 0))
    mid = side_guess // 2
    draw_helix(ImageDraw.Draw(source), mid, mid, **COMPACT_RING)
    draw_dial(ImageDraw.Draw(source), mid, mid)
    bounds = source.getchannel('A').getbbox()
elif MARK_STYLE == 'clock':
    # Crop the dial layer, so no helix dots survive in the corners.
    source = clock_layer
    bounds = clock_layer.getchannel('A').getbbox()
else:
    # text_y is the top of the wordmark, so the mark is everything above it.
    source = img
    bounds = img.crop((0, 0, width, text_y)).getchannel('A').getbbox()

if bounds is None:
    raise SystemExit(f"no artwork found for MARK_STYLE={MARK_STYLE!r}")
left, top, right, bottom = bounds
mid_x, mid_y = (left + right) / 2, (top + bottom) / 2
side = max(right - left, bottom - top) + 2 * MARK_PADDING

# Keep the square inside the canvas; the artwork is centred, so this only ever
# trims the padding rather than the mark itself.
# Clamp against the source canvas, which is not the 600x600 one for every style.
side = int(min(side, source.width, source.height))
crop_left = int(round(min(max(mid_x - side / 2, 0), source.width - side)))
crop_top = int(round(min(max(mid_y - side / 2, 0), source.height - side)))

mark = source.crop((crop_left, crop_top, crop_left + side, crop_top + side))
mark.save(MARK_PATH)
print(f"clock mark -> {MARK_PATH}  ({side}x{side} at +{crop_left}+{crop_top})")

if os.path.isdir(os.path.dirname(PYTHON_MARK_PATH)):
    shutil.copyfile(MARK_PATH, PYTHON_MARK_PATH)
    print(f"clock mark -> {PYTHON_MARK_PATH}")
else:
    print(f"skipped {PYTHON_MARK_PATH} (directory not found)")

print("\nNext: regenerate the icon set - 7 files, never edited by hand:")
print("  Rscript -e 'pkgdown::build_favicons(overwrite = TRUE)'")