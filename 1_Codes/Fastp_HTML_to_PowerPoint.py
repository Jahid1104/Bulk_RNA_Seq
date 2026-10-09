##Bulk RNA Seq Data Analysis Workflow by Md Jahid Hasan Jone##

# fastp HTML reports -> PowerPoint (charts as images, one slide per sample).
# Run after 1_Bulk_RNASeq_fastp.sh.
#
# Reads every .html file in 5_Results/1_fastp/ and writes one slide per file:
#
#     5_Results/1_fastp/fastp_charts.pptx
#
# fastp HTML reports draw their charts live in the browser with Plotly.js, so
# there are no embedded image files to extract. This script opens each HTML
# file in a headless browser (Chromium), lets the charts draw, and screenshots
# each chart. All chart images from one HTML file are placed together on one slide.
# Each slide is titled with the sample's own file name (e.g. 1F_T0_24h_R1),
# so samples from different tissues can be processed together in one run.
#
# The chart titles and positions below assume paired-end reports (11 charts).
#
# Needs (once):
#     pip install playwright python-pptx pillow
#     playwright install chromium
# On Linux, if Chromium fails to start, also run: playwright install-deps chromium

import os
import re
import sys
import shutil
import asyncio
import tempfile
from pathlib import Path
from playwright.async_api import async_playwright
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from PIL import Image


# =====================================================================
# PATHS  <- EDIT THESE (replace /.../.../ with your path)
# =====================================================================
html_dir = "/.../.../Bulk_RNA_seq/5_Results/1_fastp"
output_pptx = "/.../.../Bulk_RNA_seq/5_Results/1_fastp/fastp_charts.pptx"

KEEP_CHART_IMAGES = False   # <- EDIT: True keeps the chart PNGs in chart_images/ next to the .pptx
DEBUG = False               # <- EDIT: True prints the size and position of every chart


# CSS selector that matches every chart container in a fastp report.
# fastp draws each chart inside <div class="figure" id="plot_...">
CHART_SELECTOR = "div.figure"


async def render_file(page, html_path, out_subdir):
    os.makedirs(out_subdir, exist_ok=True)
    # as_uri() builds a valid file URL on Windows (Z:/...), Mac and Linux
    await page.goto(Path(html_path).resolve().as_uri())

    # Give Plotly time to finish drawing all charts
    await page.wait_for_timeout(1500)
    try:
        await page.wait_for_selector(CHART_SELECTOR, timeout=5000)
    except Exception:
        print(f"  ! No charts found in {html_path}, skipping.")
        return []

    charts = await page.query_selector_all(CHART_SELECTOR)
    saved = []
    for i, chart in enumerate(charts):
        out_path = os.path.join(out_subdir, f"chart_{i:02d}.png")
        try:
            await chart.screenshot(path=out_path)
            saved.append(out_path)
        except Exception as e:
            print(f"  ! Could not screenshot chart {i} in {html_path}: {e}")
    return saved


async def render_all(html_files, img_dir):
    results = {}  # html filename -> list of chart image paths
    async with async_playwright() as p:
        browser = await p.chromium.launch()
        page = await browser.new_page(viewport={"width": 1000, "height": 700})

        for html_path in html_files:
            base = os.path.splitext(os.path.basename(html_path))[0]
            out_subdir = os.path.join(img_dir, base)
            print(f"Rendering {base} ...")
            imgs = await render_file(page, html_path, out_subdir)
            results[base] = imgs
            print(f"  -> {len(imgs)} chart(s) captured")

        await browser.close()
    return results


# =====================================================================
# SLIDE LAYOUT  <- EDIT ANYTHING HERE
# =====================================================================
prs = Presentation()
prs.slide_width = Inches(13.333)   # 16:9 widescreen
prs.slide_height = Inches(7.5)
blank_layout = prs.slide_layouts[6]  # fully blank layout

SLIDE_W = prs.slide_width
SLIDE_H = prs.slide_height
MARGIN = Inches(0.3)
TITLE_H = Inches(0.5)

# Define figure titles globally (0-indexed, corresponding to chart_00.png, chart_01.png, etc.)
ALL_FIGURE_TITLES = [
    "Insert Size Distribution",
    "Read1 Quality Before Filtering",
    "Read1 Quality After Filtering",
    "Read2 Quality Before Filtering",
    "Read2 Quality After Filtering",
    "Quality Score Before Filtering",
    "Quality Score After Filtering",
    "Read1 Base Contents Before Filtering",
    "Read1 Base Contents After Filtering",
    "Read2 Base Contents Before Filtering",
    "Read2 Base Contents After Filtering"
]

# Map 0-indexed chart number to (target_row, target_col) in 0-indexed form
CHART_LAYOUT_MAP = {
    0: (0, 0),  # Figure 1: Insert Size Distribution
    5: (0, 1),  # Figure 6: Quality Score Before Filtering
    6: (0, 2),  # Figure 7: Quality Score After Filtering
    1: (1, 0),  # Figure 2: Read1 Quality Before Filtering
    2: (1, 1),  # Figure 3: Read1 Quality After Filtering
    3: (1, 2),  # Figure 4: Read2 Quality Before Filtering
    4: (1, 3),  # Figure 5: Read2 Quality After Filtering
    7: (2, 0),  # Figure 8: Read1 Base Contents Before Filtering
    8: (2, 1),  # Figure 9: Read1 Base Contents After Filtering
    9: (2, 2),  # Figure 10: Read2 Base Contents Before Filtering
    10: (2, 3)  # Figure 11: Read2 Base Contents After Filtering
}

# Explicit chart positions in inches (top-left corner)
CHART_POSITION_MAP_INCHES = {
    (0, 0): (0.5, 0.5),
    (0, 1): (4.68, 0.5),
    (0, 2): (8.85, 0.5),
    (1, 0): (0.5, 2.69),
    (1, 1): (3.69, 2.69),
    (1, 2): (6.87, 2.69),
    (1, 3): (10.05, 2.69),
    (2, 0): (0.5, 5.0),
    (2, 1): (3.69, 5.0),
    (2, 2): (6.87, 5.0),
    (2, 3): (10.05, 5.0)
}

# Estimated height for figure titles (adjust as needed)
FIGURE_TITLE_HEIGHT = Inches(0.2)
INTERNAL_PADDING = Inches(0.1)  # Padding within the calculated cell for the image itself


def get_next_coord(current_r, current_c, current_coord_val, is_horizontal, position_map, slide_max_dim_in_inches, margin_in_inches):
    """
    Calculates the extent (width or height) of a cell based on explicit start positions.
    Args:
        current_r (int): Current row index.
        current_c (int): Current column index.
        current_coord_val (float): The x (if is_horizontal) or y (if not is_horizontal) start of the current cell.
        is_horizontal (bool): True for width (x-coordinates), False for height (y-coordinates).
        position_map (dict): The CHART_POSITION_MAP_INCHES.
        slide_max_dim_in_inches (float): Total slide width or height in inches.
        margin_in_inches (float): The margin to consider from the slide edge.
    Returns:
        float: The calculated width or height of the cell.
    """
    if is_horizontal:  # Calculating width, looking for next X
        # Try to find the next column in the same row
        for (r, c), (x, _) in position_map.items():
            if r == current_r and c == current_c + 1:
                return x - current_coord_val
        # If no next column in same row, it's the last column, extend to slide end
        # (Slide width - current x - right margin)
        return slide_max_dim_in_inches - current_coord_val - margin_in_inches
    else:  # Calculating height, looking for next Y
        # Try to find the first column of the next row
        for (r, c), (_, y) in position_map.items():
            if r == current_r + 1 and c == 0:  # Assuming first col of next row defines its y-start
                return y - current_coord_val
        # If no next row, it's the last row, extend to slide end
        # (Slide height - current y - bottom margin)
        return slide_max_dim_in_inches - current_coord_val - margin_in_inches


def add_slide_with_grid(slide_title, image_paths, all_figure_titles):
    slide = prs.slides.add_slide(blank_layout)

    # Slide title
    title_box = slide.shapes.add_textbox(MARGIN, Inches(0.05), SLIDE_W - 2 * MARGIN, TITLE_H)
    tf = title_box.text_frame
    tf.text = slide_title
    tf.paragraphs[0].font.size = Pt(20)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.name = 'Times New Roman'

    if not image_paths:
        return

    charts_to_place = []
    for idx_in_list, img_path in enumerate(image_paths):
        try:
            original_chart_idx = int(os.path.basename(img_path).split('_')[1].split('.')[0])
        except (ValueError, IndexError):
            print(f"Warning: Could not parse chart index from {img_path}. Skipping.")
            continue

        if original_chart_idx in CHART_LAYOUT_MAP and original_chart_idx < len(all_figure_titles):
            charts_to_place.append((original_chart_idx, img_path, all_figure_titles[original_chart_idx]))

    charts_to_place.sort(key=lambda x: CHART_LAYOUT_MAP[x[0]])

    for original_chart_idx, img_path, fig_title_text in charts_to_place:
        r, c = CHART_LAYOUT_MAP[original_chart_idx]

        # Verify image file existence and size
        if not os.path.exists(img_path):
            print(f"ERROR: Image file not found at {img_path}")
            continue  # Skip to next image

        if os.path.getsize(img_path) == 0:
            print(f"ERROR: Image file at {img_path} is empty (0 bytes).")
            continue  # Skip to next image

        start_x_inches, start_y_inches = CHART_POSITION_MAP_INCHES[(r, c)]

        # Calculate max_cell_w_inches
        max_cell_w_inches = get_next_coord(
            r, c, start_x_inches, True, CHART_POSITION_MAP_INCHES, SLIDE_W.inches, MARGIN.inches
        )

        # Calculate max_cell_h_inches
        max_cell_h_inches = get_next_coord(
            r, c, start_y_inches, False, CHART_POSITION_MAP_INCHES, SLIDE_H.inches, MARGIN.inches
        )

        # Apply fixed dimensions based on row
        if r == 0:  # First row
            final_w_inches = 4.0
            final_h_inches = 2.0
        elif r == 1 or r == 2:  # Second and third rows
            final_w_inches = 3.0
            final_h_inches = 2.0
        else:
            # Fallback or error handling for unexpected rows
            print(f"Warning: Unexpected row index {r}. Using default scaling.")
            with Image.open(img_path) as im:
                iw, ih = im.size  # Original image dimensions in pixels

            # Subtract space for internal padding and title height from the effective image dimensions
            effective_max_w_for_image = max_cell_w_inches - 2 * INTERNAL_PADDING.inches
            effective_max_h_for_image = max_cell_h_inches - 2 * INTERNAL_PADDING.inches - FIGURE_TITLE_HEIGHT.inches

            # Ensure calculated max dimensions are positive
            effective_max_w_for_image = max(0.1, effective_max_w_for_image)
            effective_max_h_for_image = max(0.1, effective_max_h_for_image)

            # Calculate scale to fit within the effective max dimensions (ratio of inches/pixel)
            scale_to_fit_width = effective_max_w_for_image / iw
            scale_to_fit_height = effective_max_h_for_image / ih
            optimal_scale = min(scale_to_fit_width, scale_to_fit_height)

            # Calculate final image dimensions in inches
            final_w_inches = iw * optimal_scale
            final_h_inches = ih * optimal_scale

        # Ensure dimensions are at least a minimal value to prevent rendering issues
        w = Inches(max(0.01, final_w_inches))
        h = Inches(max(0.01, final_h_inches))

        if DEBUG:
            print(f"DEBUG: Chart {original_chart_idx:02d} on slide '{slide_title}': {w.inches:.3f}x{h.inches:.3f} inches at position ({start_x_inches:.3f}, {start_y_inches:.3f})")

        # Center the image within its allocated cell space defined by max_cell_w_inches and max_cell_h_inches
        img_final_x = Inches(start_x_inches + (max_cell_w_inches - w.inches) / 2)
        img_final_y = Inches(start_y_inches + (max_cell_h_inches - h.inches - FIGURE_TITLE_HEIGHT.inches) / 2)

        slide.shapes.add_picture(img_path, img_final_x, img_final_y, width=w, height=h)

        # Add figure title below the image, centered
        # The title box width should be the image width for centering relative to image
        title_box_x = img_final_x
        title_box_y = img_final_y + h + INTERNAL_PADDING.inches  # Small gap between image and title
        title_box_w = w

        fig_title_shape = slide.shapes.add_textbox(
            title_box_x, title_box_y, title_box_w, FIGURE_TITLE_HEIGHT
        )
        tf_fig = fig_title_shape.text_frame
        tf_fig.text = fig_title_text
        tf_fig.paragraphs[0].font.size = Pt(9)
        tf_fig.paragraphs[0].alignment = PP_ALIGN.CENTER
        tf_fig.paragraphs[0].font.name = 'Times New Roman'
        tf_fig.word_wrap = False


# Sort slides in natural order (1F, 2F, ..., 10F, 1L, 2L, ...) instead of
# plain alphabetical order (which would put "10F" before "2F").
def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def main():
    if "..." in html_dir:
        sys.exit("Edit the paths at the top of this script (replace /.../.../ with your path).")

    html_files = sorted(
        (os.path.join(html_dir, f) for f in os.listdir(html_dir) if f.lower().endswith(".html")),
        key=natural_key,
    )
    if not html_files:
        sys.exit(f"No .html files found in {html_dir}")
    print(f"Found {len(html_files)} HTML file(s).")

    if KEEP_CHART_IMAGES:
        img_dir = os.path.join(os.path.dirname(output_pptx), "chart_images")
        os.makedirs(img_dir, exist_ok=True)
    else:
        img_dir = tempfile.mkdtemp(prefix="fastp_chart_images_")

    try:
        chart_results = asyncio.run(render_all(html_files, img_dir))

        # Each HTML file already carries its own sample name (tissue, temp,
        # time, replicate, etc.), so it is used directly as the slide title.
        for base in sorted(chart_results.keys(), key=natural_key):
            imgs = chart_results[base]
            add_slide_with_grid(base, sorted(imgs), ALL_FIGURE_TITLES)

        os.makedirs(os.path.dirname(output_pptx), exist_ok=True)
        prs.save(output_pptx)
        print(f"\nDone. Saved {len(chart_results)} slide(s) to {output_pptx}")
    finally:
        if not KEEP_CHART_IMAGES:
            shutil.rmtree(img_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
