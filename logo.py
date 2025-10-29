from PIL import Image, ImageDraw, ImageFont
import numpy as np

# Create image
width, height = 800, 800
img = Image.new('RGB', (width, height), color='white')
draw = ImageDraw.Draw(img)

# Color scheme - scientific blue gradient
bg_color = (245, 250, 255)
primary_color = (30, 100, 180)
accent_color = (100, 180, 230)

# Fill background
img = Image.new('RGB', (width, height), color=bg_color)
draw = ImageDraw.Draw(img)

# Draw circular gradient background
center_x, center_y = width // 2, height // 2
for r in range(350, 0, -5):
    alpha = r / 350
    color = tuple(int(bg_color[i] + (primary_color[i] - bg_color[i]) * (1 - alpha)) for i in range(3))
    draw.ellipse([center_x - r, center_y - r, center_x + r, center_y + r], fill=color)

# Draw DNA helix-inspired spiral (adjusted radius for better centering and less overlap)
num_points = 100
for i in range(num_points):
    angle = (i / num_points) * 4 * np.pi
    radius = 150 + i * 0.6  # Start larger and spiral out more controllably
    x1 = center_x + radius * np.cos(angle)
    y1 = center_y + radius * np.sin(angle)
    x2 = center_x + radius * np.cos(angle + np.pi)
    y2 = center_y + radius * np.sin(angle + np.pi)
    
    size = int(4 + i * 0.08)  # Slightly smaller points to reduce clutter
    color_val = int(100 + (i / num_points) * 130)
    point_color = (color_val, color_val + 50, 230)
    
    draw.ellipse([x1-size, y1-size, x1+size, y1+size], fill=point_color)
    draw.ellipse([x2-size, y2-size, x2+size, y2+size], fill=accent_color)

# Draw clock elements - representing "age" (centered for better balance)
clock_radius = 120
clock_center = (center_x, center_y)  # Centered now for proper placement

# Clock circle
draw.ellipse([clock_center[0] - clock_radius, clock_center[1] - clock_radius,
              clock_center[0] + clock_radius, clock_center[1] + clock_radius],
             outline=primary_color, width=8)

# Clock hands (adjusted angles slightly for a more balanced "time" representation)
hand_angle = -np.pi / 4  # Adjusted to about 7:30 position for better visual flow
hand_length = 80
hand_x = clock_center[0] + hand_length * np.cos(hand_angle)
hand_y = clock_center[1] + hand_length * np.sin(hand_angle)
draw.line([clock_center[0], clock_center[1], hand_x, hand_y], fill=primary_color, width=6)

minute_angle = np.pi / 4  # Adjusted to complement the hour hand
minute_length = 100
minute_x = clock_center[0] + minute_length * np.cos(minute_angle)
minute_y = clock_center[1] + minute_length * np.sin(minute_angle)
draw.line([clock_center[0], clock_center[1], minute_x, minute_y], fill=accent_color, width=4)

# Center dot
draw.ellipse([clock_center[0]-10, clock_center[1]-10, 
              clock_center[0]+10, clock_center[1]+10], fill=primary_color)

# Add "tAge" text (adjusted positioning for better vertical balance)
try:
    # Try to use a nice font
    font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 100)
    font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 40)
except:
    font_large = ImageFont.load_default()
    font_small = ImageFont.load_default()

# Main text
text = "tAge"
bbox = draw.textbbox((0, 0), text, font=font_large)
text_width = bbox[2] - bbox[0]
text_x = (width - text_width) // 2
text_y = height - 180  # Slightly higher for better centering with centered clock

# Add shadow
shadow_offset = 3
draw.text((text_x + shadow_offset, text_y + shadow_offset), text, 
          fill=(200, 200, 200), font=font_large)
draw.text((text_x, text_y), text, fill=primary_color, font=font_large)

# Subtitle
subtitle = "Transcriptomic Age"
bbox_sub = draw.textbbox((0, 0), subtitle, font=font_small)
sub_width = bbox_sub[2] - bbox_sub[0]
sub_x = (width - sub_width) // 2
sub_y = text_y + 120  # Adjusted to fit better

draw.text((sub_x, sub_y), subtitle, fill=accent_color, font=font_small)

# Save
img.save('tAge_logo.png')
print("Logo saved to tAge_logo.png")