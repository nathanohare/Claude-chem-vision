"""Mock substrate-scope table PDF: vector structures (RDKit SVG) with bold
compound IDs and yield/ee/dr text underneath, like a journal scope figure.
Page 1 has a text layer; page 2 is the same table rasterised (no text layer)."""
import json, io, os
os.makedirs("images", exist_ok=True)
import pymupdf as fitz
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
rdDepictor.SetPreferCoordGen(True)

# Asymmetric ketone reduction scope (illustrative values).
SCOPE = [
 ("3a", "O[C@@H](C)c1ccccc1", ["92%, 95% ee"]),
 ("3b", "O[C@@H](C)c1ccc(OC)cc1", ["88%, 97:3 er"]),
 ("3c", "O[C@@H](C)c1ccc(Br)cc1", ["90% yield", "94% ee"]),
 ("3d", "O[C@@H](CC)c1ccccc1", ["71%, 89% ee", "(48 h)"]),
 ("3e", "O[C@@H](C)c1ccc2ccccc2c1", ["85%, >99% ee"]),
 ("3f", "O[C@H]1CCCc2ccccc21", ["64%, 91% ee", ">20:1 dr"]),
]

def svg(smi):
    d = rdMolDraw2D.MolDraw2DSVG(220, 170)
    d.drawOptions().useBWAtomPalette()
    rdMolDraw2D.PrepareAndDrawMolecule(d, Chem.MolFromSmiles(smi))
    d.FinishDrawing()
    return d.GetDrawingText()

doc = fitz.open()
page = doc.new_page(width=612, height=792)
page.insert_text((54, 60), "Table 2. Substrate scope of the asymmetric transfer hydrogenation", fontsize=11, fontname="helv")
page.insert_text((54, 76), "Conditions: 1 (0.5 mmol), Ru catalyst (1 mol%), HCO2H/Et3N, 30 °C, 24 h. Isolated yields; ee by HPLC.", fontsize=8, fontname="helv")
boxes = {}
for i, (cid, smi, lines) in enumerate(SCOPE):
    col, row = i % 3, i // 3
    x, y = 54 + col * 170, 100 + row * 230
    s = fitz.open("svg", svg(smi).encode())
    pdfbytes = s.convert_to_pdf(); s.close()
    src = fitz.open("pdf", pdfbytes)
    r = fitz.Rect(x, y, x + 160, y + 124)
    page.show_pdf_page(r, src, 0)
    boxes[cid] = list(r)
    page.insert_text((x + 60, y + 140), cid, fontsize=10, fontname="hebo")
    for j, t in enumerate(lines):
        page.insert_text((x + 40, y + 154 + 12 * j), t, fontsize=9, fontname="helv")
# Page 2: rasterised copy (no text layer), like a scanned/bitmap figure.
pix = page.get_pixmap(dpi=200)
p2 = doc.new_page(width=612, height=792)
p2.insert_image(p2.rect, stream=pix.tobytes("png"))
doc.save("images/scope_table.pdf")
json.dump({cid: {"smiles": smi, "labels": lines, "box_pt": boxes[cid]} for cid, smi, lines in SCOPE},
          open("scope_truth.json", "w"), indent=1)
print("ok")
