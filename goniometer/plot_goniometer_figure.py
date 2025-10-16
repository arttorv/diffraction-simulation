import os
from PIL import Image
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Run this file to plot simulation/observed images comparison from goniometer. 
# The experimental images are taken from "P600_images" and "P2100_images"
# Calculated images are taken from the "export"-folder and are created by program 2.

# Please insert values below: 

# ------------------ #
patch = 'P600' # P2100 or P600
simulation_resolution = 0.1 # Degrees (0.1 - 1)
# ------------------ # 

if patch == 'P2100':
    folders = {
        '0%': 'goniometer/P2100_images/extracted_ROIs_0',
        '5%': 'goniometer/P2100_images/extracted_ROIs_5',
        '10%': 'goniometer/P2100_images/extracted_ROIs_10'
    }
    # Periods of strained grating 
    sim_p = {
        '0%': '2150',
        '5%': '2258',
        '10%': '2365'
    }
    frame_stop = 10 # Limit to observed angle in experiment

elif patch == "P600":
    folders = {
        '0%': 'goniometer/P600_images/extracted_rois_0',
        '5%': 'goniometer/P600_images/extracted_rois_5',
        '10%': 'goniometer/P600_images/extracted_rois_10'
    }
    # Periods of strained grating 
    sim_p = {
        '0%': '630',
        '5%': '660',
        '10%': '690'
    }
    frame_stop = 20 # Limit to observed angle in experiment

# Common frame numbers and folders for each strain
frame_numbers = list(range(90, frame_stop, -1))  # From 90 to 10

# Load and rotate images per strain
all_images = {}
max_widths = []
max_heights = []

for strain_label, folder in folders.items():
    images = []
    for i in frame_numbers:
        path = os.path.join(folder, f'frame_{i}.png')
        img = Image.open(path).rotate(90, expand=True)
        # crop the images to square size from the center
        min_size = min(img.size)
        # crop the first 50 images from center, the rest from bottom
        if i < 37:
            # crop from bottom
            left = (img.width - min_size) / 2
            top = img.height - min_size
            right = (img.width + min_size) / 2
            bottom = img.height *2
            # img = img.crop((left, top, right, bottom))
        else:
            left = (img.width - min_size) / 2
            top = (img.height - min_size) / 2
            right = (img.width + min_size) / 2
            bottom = (img.height + min_size) / 2
            # img = img.crop((left, top, right, bottom))
        images.append(img)
    all_images[strain_label] = images
    widths, heights = zip(*(img.size for img in images))
    max_widths.append(sum(widths))
    max_heights.append(max(heights))

# Determine canvas size
total_width = max(max_widths)
image_height = max(max_heights)
spacing = 3*image_height  # vertical spacing between strains
n_strains = len(folders)
total_height = n_strains * image_height + (n_strains - 1) * spacing + 2*image_height  # + image_height for last simulation bar

# Create canvas
fig, ax = plt.subplots(figsize=(total_width / 300, total_height / 100))
ax.axis('off')

# Plot each strain's images with vertical offset
strain_order = ['0%', '5%', '10%']  # top to bottom
for idx, strain in enumerate(strain_order):
    images = all_images[strain]
    y_offset = total_height - (idx + 1) * image_height - idx * spacing
    x_offset = 0

    for i, img in zip(frame_numbers, images):
        ax.imshow(img, extent=[x_offset, x_offset + img.width, y_offset, y_offset + img.height])
        x_offset += img.width

    # Add strain label to left
    if idx == 0:
        ax.text(-100, y_offset/2 + image_height / 2, "Strain", ha='right', va='center', fontsize=14, rotation=90)
    ax.text(-30, y_offset + image_height / 2, strain, ha='right', va='center', fontsize=12, rotation=90)

    # Only add angle labels to the bottom row (10%)
    if strain == '10%':
        x_offset = 0
        for i in frame_numbers:
            img = all_images[strain][frame_numbers.index(i)]
            center_x = x_offset + img.width / 2
            angle = str(90 + 4 - i)  # adjust offset if needed
            # add the angle label only for 10, 20, 30, ..., 90
            if (i-4) % 10 == 0:
                ax.text(center_x, y_offset - 200, angle, ha='center', va='center', fontsize=12)
            else:
                ax.text(center_x, y_offset - 120, '', ha='center', va='center', fontsize=6)
            x_offset += img.width

# SIMULATION COLOR BARS
# ---------------------
angle_res_scale = simulation_resolution
simulation_periods = [sim_p['0%'], sim_p['5%'], sim_p['10%']]
img_width = all_images['10%'][0].size[0] * angle_res_scale*1.01
bar_height = image_height/2
sim_y_offset = image_height +20 # Adjust vertical offset to place under image row
angle_frames = [round(a, 1) for a in np.arange(90, frame_stop - angle_res_scale, -angle_res_scale)]

x_coords = [(i * img_width, (i + 1) * img_width) for i in range(len(angle_frames))]
# ---------------------

# Loop over your simulation CSV files
count = 1

for idx, period in enumerate(simulation_periods):
    # Load color data: space-separated, no headers
    filepath = f'export/diffraction_data_{period}.csv'
    sim_df = pd.read_csv(filepath, sep=",", header=None, names=["angle", "r", "g", "b"])

    # Convert angle to int, RGB to float in [0, 1]
    sim_df["angle"] = sim_df["angle"].astype(float)
    sim_df[["r", "g", "b"]] = sim_df[["r", "g", "b"]].astype(float)
    color_data = {row.angle: (row.r, row.g, row.b) for _, row in sim_df.iterrows()}

    # Plot color squares aligned to bottom images
    y_base = total_height - (idx + 1) * image_height - idx * spacing - sim_y_offset
    for img_idx, frame_num in enumerate(angle_frames):
        angle = 90 + 4 - frame_num
        angle = round(angle, 1)  # round to match the CSV
        color = color_data.get(angle, (0, 0, 0))  # default black
        x0, x1 = x_coords[img_idx]
        ax.add_patch(plt.Rectangle((x0, y_base), x1 - x0, bar_height, color=color, linewidth=0))
    

    ax.text(total_width + 110, y_base + bar_height / 2, f'Sim', ha='right', va='center', fontsize=12)
    ax.text(total_width + 110, y_base + 1.8 * image_height , f'Exp', ha='right', va='center', fontsize=12)
    count += 1 

# Add global label below all
# ax.text(total_width / 2, y_base - 20, "Simulated diffraction color (per period)", ha='center', va='center', fontsize=9)


# Final layout
ax.text(total_width/2, -130, "Angle (degrees)", ha='center', va='center', fontsize=12)

# ax.text(total_width/2, total_height + 50, "Diffracted color from ", ha='center', va='center', fontsize=12)
ax.set_xlim(0, total_width)
ax.set_ylim(0, total_height)
plt.tight_layout()
plt.show()
