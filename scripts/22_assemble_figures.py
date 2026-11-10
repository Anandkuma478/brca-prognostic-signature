"""Assemble the multi-panel manuscript figures from the single panels written by the R scripts.

Figure 9  (fig_km_panels.png):              km_*.png from 12_survival_evaluation.R, 3 columns
Figure 10 (fig_calibration_panels.png):     calibration_*.png from 12_survival_evaluation.R, 3 columns
Figure 11 (fig_batch_pca_before_after.png): batch_panel_*.png from 15_batch_correction_sensitivity.R, stacked

Each figure is resized to 650 pixels per inch of its printed width in the manuscript, so it is above
600 dpi when placed. The panels themselves are drawn at 600 dpi, so no panel is enlarged.
Run from the project root, after steps 13 and 15:  python3 scripts/22_assemble_figures.py
"""
import os

from PIL import Image

F = "figures/"
PPI = 650


def grid(names, cols, out, print_width_in):
    ims = [Image.open(F + n).convert("RGB") for n in names]
    cw = max(i.width for i in ims)
    ch = max(i.height for i in ims)
    rows = (len(ims) + cols - 1) // cols
    sheet = Image.new("RGB", (cols * cw, rows * ch), "white")
    for k, im in enumerate(ims):
        sheet.paste(im, ((k % cols) * cw, (k // cols) * ch))
    w = round(print_width_in * PPI)
    assert w <= sheet.width, f"{out}: panels too small to reach {PPI} ppi without enlarging"
    sheet = sheet.resize((w, round(sheet.height * w / sheet.width)), Image.LANCZOS)
    sheet.save(F + out, optimize=True, dpi=(PPI, PPI))
    print(out, sheet.size, f"{os.path.getsize(F + out) / 1e6:.1f} MB")


grid(["km_SCANB_modelA.png", "km_TCGABRCA_modelA.png", "km_METABRIC_modelA.png",
      "km_SCANB_modelB.png", "km_TCGABRCA_modelB.png", "km_METABRIC_modelB.png"], 3, "fig_km_panels.png", 6.0)
grid(["calibration_SCANB_modelA.png", "calibration_TCGABRCA_modelA.png", "calibration_METABRIC_modelA.png",
      "calibration_SCANB_modelB.png", "calibration_TCGABRCA_modelB.png", "calibration_METABRIC_modelB.png"],
     3, "fig_calibration_panels.png", 6.0)
grid(["batch_panel_scanb.png", "batch_panel_tcga.png", "batch_panel_metabric.png"], 1, "fig_batch_pca_before_after.png", 5.6)
