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
for i in range(num_points):
    angle = (i / num_points) * 3 * np.pi
    radius = 120 + i * 1.2
    x1 = center_x + radius * np.cos(angle)
    y1 = center_y + radius * np.sin(angle) - 50
    x2 = center_x + radius * np.cos(angle + np.pi)
    y2 = center_y + radius * np.sin(angle + np.pi) - 50
    
    size = int(5 + i * 0.06)
    alpha = int(180 + (i / num_points) * 75)
    
    draw.ellipse([x1-size, y1-size, x1+size, y1+size], 
                 fill=primary_color + (alpha,))
    draw.ellipse([x2-size, y2-size, x2+size, y2+size], 
                 fill=accent_color + (alpha,))

# Simplified clock - cleaner design
clock_radius = 100
clock_y = center_y - 50

# Clock circle
draw.ellipse([center_x - clock_radius, clock_y - clock_radius,
              center_x + clock_radius, clock_y + clock_radius],
             outline=primary_color + (255,), width=6)

# Clock hands
hand_angle = -np.pi / 3
hand_length = 60
hand_x = center_x + hand_length * np.cos(hand_angle)
hand_y = clock_y + hand_length * np.sin(hand_angle)
draw.line([center_x, clock_y, hand_x, hand_y], 
          fill=primary_color + (255,), width=5)

minute_angle = np.pi / 6
minute_length = 80
minute_x = center_x + minute_length * np.cos(minute_angle)
minute_y = clock_y + minute_length * np.sin(minute_angle)
draw.line([center_x, clock_y, minute_x, minute_y], 
          fill=accent_color + (255,), width=4)

# Center dot
draw.ellipse([center_x-8, clock_y-8, center_x+8, clock_y+8], 
             fill=primary_color + (255,))

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

# Save as PNG with transparency
img.save('tAge_logo.png')
print("Transparent logo saved to tAge_logo.png")