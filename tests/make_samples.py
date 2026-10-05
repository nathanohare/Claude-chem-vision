"""Builds the synthetic test set: varied drawing styles, degraded scans,
a reaction scheme, and a mock journal page PDF. Ground truth in truth.json."""
import io, json, os, random
os.makedirs("images", exist_ok=True)
from PIL import Image, ImageDraw, ImageFilter, ImageFont
from rdkit import Chem
from rdkit.Chem import rdAbbreviations
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)  # closer to ChemDraw-style layouts for cages/macrocycles

random.seed(0)
TRUTH = {
 "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
 "testosterone": "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O",
 "atorvastatin": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
 "boc_ester": "COC(=O)[C@H](Cc1ccc(OC)cc1)NC(=O)OC(C)(C)C",
 "paclitaxel": "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C",
 "strychnine": "O=C7N2c1ccccc1[C@@]64[C@@H]2[C@@H]3[C@@H](OC/C=C5\\[C@@H]3C[C@@H]6N(CC4)C5)C7",
 "caffeine_rot": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
 "sildenafil_scan": "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1",
}

def draw(smi, size=(500, 400), bw=False, width=2, abbrev=False, font=0.0):
    m = Chem.MolFromSmiles(smi)
    if abbrev:
        m = rdAbbreviations.CondenseMolAbbreviations(m, rdAbbreviations.GetDefaultAbbreviations())
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    o = d.drawOptions()
    if bw: o.useBWAtomPalette()
    o.bondLineWidth = width
    if font: o.minFontSize = int(font)
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText())).convert("RGB")

def jpeg(im, q):
    b = io.BytesIO(); im.save(b, "JPEG", quality=q); return Image.open(b).convert("RGB")

def scanify(im):
    g = im.convert("L").filter(ImageFilter.GaussianBlur(0.8))
    px = g.load()
    for _ in range(g.width * g.height // 40):
        x, y = random.randrange(g.width), random.randrange(g.height)
        px[x, y] = random.choice([0, 120, 255])
    return jpeg(g.rotate(1.5, fillcolor=255, expand=True).convert("RGB"), 35)

draw(TRUTH["aspirin"]).save("images/aspirin.png")
draw(TRUTH["testosterone"], bw=True, width=3).save("images/testosterone.png")
jpeg(draw(TRUTH["atorvastatin"], size=(320, 260), bw=True), 40).save("images/atorvastatin.jpg")
draw(TRUTH["boc_ester"], bw=True, abbrev=True).save("images/boc_ester.png")
draw(TRUTH["paclitaxel"], size=(800, 650), bw=True).save("images/paclitaxel.png")
draw(TRUTH["strychnine"], size=(500, 450), bw=True).save("images/strychnine.png")
draw(TRUTH["caffeine_rot"], bw=True).rotate(25, fillcolor="white", expand=True).save("images/caffeine_rot.png")
scanify(draw(TRUTH["sildenafil_scan"], size=(600, 450), bw=True)).save("images/sildenafil_scan.jpg")

# Reaction scheme: two-step sequence drawn the way papers do (reagents over arrows).
font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 17)
def arrow(w, top, bottom):
    im = Image.new("RGB", (w, 380), "white"); d = ImageDraw.Draw(im)
    d.line([(10, 190), (w - 12, 190)], fill="black", width=3)
    d.polygon([(w - 4, 190), (w - 20, 181), (w - 20, 199)], fill="black")
    for i, t in enumerate(top): d.text((14, 140 - 22 * (len(top) - 1 - i)), t, fill="black", font=font)
    for i, t in enumerate(bottom): d.text((14, 205 + 22 * i), t, fill="black", font=font)
    return im
def label(im, t):
    d = ImageDraw.Draw(im); d.text((im.width // 2 - 8, im.height - 40), t, fill="black", font=font); return im
SCHEME = {
 "1": "COc1ccc(Br)cc1", "2": "OB(O)c1ccc(C=O)cc1", "3": "COc1ccc(-c2ccc(C=O)cc2)cc1", "4": "COc1ccc(-c2ccc(CO)cc2)cc1",
}
tiles = [label(draw(SCHEME["1"], (260, 380), bw=True), "1"),
         label(Image.new("RGB", (60, 380), "white"), ""),
         label(draw(SCHEME["2"], (280, 380), bw=True), "2"),
         arrow(230, ["Pd(PPh3)4 (5 mol%)", "K2CO3"], ["dioxane/H2O", "90 °C, 12 h, 87%"]),
         label(draw(SCHEME["3"], (380, 380), bw=True), "3"),
         arrow(170, ["NaBH4"], ["MeOH, 0 °C", "95%"]),
         label(draw(SCHEME["4"], (380, 380), bw=True), "4")]
d = ImageDraw.Draw(tiles[1]); d.text((20, 175), "+", fill="black", font=ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 34))
W = sum(t.width for t in tiles); sch = Image.new("RGB", (W, 380), "white"); x = 0
for t in tiles: sch.paste(t, (x, 0)); x += t.width
sch.save("images/scheme_suzuki.png")

# Mock journal page: text, a numbered compound figure, and the scheme.
import fitz
doc = fitz.open(); page = doc.new_page(width=612, height=792)
page.insert_text((54, 60), "Synthesis of biaryl alcohol 4", fontsize=15)
body = ("Aryl bromide 1 was coupled with boronic acid 2 under standard Suzuki-Miyaura conditions "
        "to give aldehyde 3, which was reduced with sodium borohydride to afford 4 (Scheme 1). "
        "Compounds 5 and 6 (Figure 1) were prepared as described previously.")
page.insert_textbox(fitz.Rect(54, 75, 558, 140), body, fontsize=10)
page.insert_text((54, 160), "Scheme 1. Two-step synthesis of 4.", fontsize=9)
b = io.BytesIO(); sch.save(b, "PNG"); page.insert_image(fitz.Rect(54, 168, 558, 168 + 504 * 380 / W), stream=b.getvalue())
fig = Image.new("RGB", (900, 420), "white")
fig.paste(label(draw(TRUTH["boc_ester"], (450, 420), bw=True, abbrev=True), "5"), (0, 0))
fig.paste(label(draw(TRUTH["testosterone"], (450, 420), bw=True), "6"), (450, 0))
b = io.BytesIO(); fig.save(b, "PNG"); page.insert_image(fitz.Rect(54, 260, 558, 260 + 504 * 420 / 900), stream=b.getvalue())
page.insert_text((54, 505), "Figure 1. Structures of 5 and 6.", fontsize=9)
doc.save("images/mock_paper.pdf")

TRUTH.update({"scheme_" + k: v for k, v in SCHEME.items()})
TRUTH["scheme_rxn"] = "COc1ccc(Br)cc1.OB(O)c1ccc(C=O)cc1>>COc1ccc(-c2ccc(C=O)cc2)cc1"
json.dump(TRUTH, open("truth.json", "w"), indent=1)
print("ok")
