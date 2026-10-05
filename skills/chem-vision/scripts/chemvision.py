#!/usr/bin/env python3
"""chemvision: helper CLI for reading chemical structures and reaction schemes
from images and PDFs.

Claude does the visual interpretation; this tool supplies what vision alone is
bad at: page rendering and zooming, a second opinion from an optical chemical
structure recognition (OCSR) engine, chemistry validation with RDKit, and
re-rendering a proposed SMILES so it can be compared against the original.

Every subcommand prints JSON to stdout so results are easy to read back.

Subcommands
  pages      Render PDF pages to PNG (and list embedded images); --grid adds overlays.
  grid       Draw a labeled coordinate grid on any image (for choosing crop boxes).
  docx-figures  Extract embedded figures (EMF/WMF/PNG/JPEG) from a .docx as PNGs.
  cdxml      Read a ChemDraw .cdxml: label text tied to each structure (+ SMILES with RDKit).
  crop       Crop (and upscale) a region of an image to look at it closely.
  recognize  Run OCSR engines (OSRA, MolScribe, DECIMER) on an image.
  build      Build a molecule from drawn atom positions and wedge/hash bonds.
  check      Validate SMILES and report formula, mass, InChIKey, groups, stereo.
  render     Draw SMILES to PNG, optionally side by side with the source image.
  compare    Score agreement between candidate SMILES (InChIKey + Tanimoto).
  reaction   Parse a reaction SMILES, check balance, and render the scheme.
  annotate   Compound IDs, yields, ee/er, dr printed next to structures.
  export     Write compound records (ID, SMILES, yield, ee, ...) to CSV/SDF/JSON.
  lookup     Look a structure up on PubChem by InChIKey or name (needs network).
  engines    Report which engines and dependencies are available.

Commands that do not need RDKit: pages, grid, crop, docx-figures, cdxml
(labels only), annotate (--text / --labels / --pdf), lookup, engines.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, DataStructs, Draw, rdFingerprintGenerator, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D

    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT, RDKIT_ERROR = True, None
except Exception as _e:  # not installed, or install blocked (e.g. network allowlist)
    HAVE_RDKIT, RDKIT_ERROR = False, f"{type(_e).__name__}: {_e}"

# Subcommands that cannot run at all without RDKit.
RDKIT_COMMANDS = {"recognize", "check", "render", "build", "compare", "reaction", "export"}
RDKIT_MISSING = {
    "error": "RDKit is not available in this environment",
    "detail": None,
    "unavailable": ["recognize (OCSR second opinion)", "check (SMILES validation, formula, MW, InChIKey)",
                    "render (redraw-and-compare)", "build (stereo from drawing)", "compare", "reaction",
                    "export", "cdxml structures (labels still work)"],
    "still_available": ["pages", "grid", "crop", "docx-figures", "cdxml (label text, values, positions)",
                        "annotate --text / --labels / --pdf", "lookup", "engines"],
    "what_to_do": "Continue by reading the structure visually and say in the report that SMILES, "
                  "stereodescriptors and formulas are unvalidated. To fix: `pip install rdkit`. If pip "
                  "reports 'No matching distribution' or PyPI returns 403, the package index is blocked by "
                  "this environment's network allowlist; retrying with other commands will not help. Run the "
                  "plugin on a machine where rdkit installs, or ask an org admin to allow PyPI.",
}


def out(obj):
    print(json.dumps(obj, indent=2, ensure_ascii=False))


# ---------------------------------------------------------------- pages / crop

def cmd_pages(a):
    try:
        import pymupdf as fitz
    except Exception:
        return pages_poppler(a)

    doc = fitz.open(a.pdf)
    os.makedirs(a.outdir, exist_ok=True)
    pages = parse_range(a.pages, len(doc)) if a.pages else range(len(doc))
    stem = os.path.splitext(os.path.basename(a.pdf))[0]
    result = []
    for i in pages:
        page = doc[i]
        pix = page.get_pixmap(dpi=a.dpi)
        path = os.path.join(a.outdir, f"{stem}_p{i + 1}.png")
        pix.save(path)
        entry = {"page": i + 1, "png": path, "size": [pix.width, pix.height]}
        if a.grid:
            entry["grid"] = draw_grid(path, path[:-4] + "_grid.png")
        if a.images:
            imgs = []
            for j, info in enumerate(page.get_images(full=True)):
                xref, w0, h0 = info[0], info[2], info[3]
                if w0 < 80 or h0 < 80:
                    continue
                try:
                    rects = page.get_image_rects(xref)
                    ip = os.path.join(a.outdir, f"{stem}_p{i + 1}_img{j + 1}.png")
                    if rects and rects[0].width > 0:
                        # Render the image as the page shows it: this applies
                        # masks and decode arrays (1-bit CCITT figures otherwise
                        # come out inverted). Native resolution, capped at 600 dpi.
                        r = rects[0]
                        dpi = min(600, max(150, int(w0 / (r.width / 72))))
                        img = page.get_pixmap(clip=r, dpi=dpi)
                    else:
                        img = fitz.Pixmap(doc, xref)
                        if img.n - img.alpha >= 4:
                            img = fitz.Pixmap(fitz.csRGB, img)
                        dpi = None
                    img.save(ip)
                    imgs.append({"png": ip, "size": [img.width, img.height], "dpi": dpi,
                                 "page_rect_pt": [round(v, 1) for v in rects[0]] if rects else None})
                except Exception as e:  # odd colourspaces etc.
                    imgs.append({"error": str(e), "xref": xref})
            entry["embedded_images"] = imgs
        if a.text:
            entry["text"] = page.get_text()
        result.append(entry)
    out({"pdf": a.pdf, "page_count": len(doc), "pages": result})


def parse_range(spec, n):
    pages = []
    for part in spec.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            pages.extend(range(int(lo) - 1, min(int(hi), n)))
        else:
            pages.append(int(part) - 1)
    return [p for p in pages if 0 <= p < n]


def cmd_crop(a):
    from PIL import Image

    im = Image.open(a.image).convert("RGB")
    W, H = im.size
    x0, y0, x1, y1 = [float(v) for v in a.box.split(",")]
    if max(x0, y0, x1, y1) <= 1.0:  # fractional coordinates
        x0, x1, y0, y1 = x0 * W, x1 * W, y0 * H, y1 * H
    pad = a.pad
    box = (max(0, int(x0) - pad), max(0, int(y0) - pad),
           min(W, int(x1) + pad), min(H, int(y1) + pad))
    c = im.crop(box)
    inverted = False
    if a.invert or (not a.no_auto_invert and is_dark(c)):
        from PIL import ImageOps
        c, inverted = ImageOps.invert(c), True
    if a.scale != 1:
        c = c.resize((int(c.width * a.scale), int(c.height * a.scale)), Image.LANCZOS)
    c.save(a.out)
    out({"crop": a.out, "box_px": list(box), "size": list(c.size), "inverted": inverted})


def is_dark(im):
    """White-on-black drawing (inverted figure extraction)?"""
    from PIL import ImageStat
    return ImageStat.Stat(im.convert("L")).mean[0] < 110


# ------------------------------------------------------------------- recognize

def engines_available():
    if HAVE_RDKIT:
        from rdkit import __version__ as rdkit_version
    else:
        rdkit_version = None
    eng = {"rdkit": rdkit_version, "osra": bool(shutil.which("osra"))}
    for name, mod in (("molscribe", "molscribe"), ("decimer", "DECIMER")):
        try:
            __import__(mod)
            eng[name] = True
        except Exception:
            eng[name] = False
    return eng


BOX_RE = re.compile(r"^(\d+)x(\d+)-(\d+)x(\d+)$")


def run_osra(image, extra=()):
    cmd = ["osra", "-f", "smi", "-p", "-b", "-c", *extra, image]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    hits = []
    for line in p.stdout.splitlines():
        parts = line.split()
        if not parts:
            continue
        h = {"engine": "osra", "smiles": parts[0]}
        nums = [x for x in parts[1:] if not BOX_RE.match(x)]
        boxes = [x for x in parts[1:] if BOX_RE.match(x)]
        if len(nums) >= 2:
            h["bond_length_px"] = float(nums[0])
            h["confidence"] = float(nums[1])
        if boxes:
            m = BOX_RE.match(boxes[0])
            h["box_px"] = [int(m.group(i)) for i in range(1, 5)]
        hits.append(h)
    return hits


def osra_ensemble(image, adaptive=False, jaggy=False, rotate=0, quick=False):
    """OSRA's output swings with preprocessing, so run several passes
    (default, adaptive threshold, thinning, median denoise, 2x upscale) and vote. Candidates
    are grouped by InChIKey; chemically invalid readings are dropped unless no
    pass produced a valid one."""
    from PIL import Image

    base = []
    if adaptive:
        base.append("-i")
    if jaggy:
        base.append("-j")
    if rotate:
        base += ["-R", str(int(rotate))]
    tmp = tempfile.mkdtemp()
    im0 = Image.open(image).convert("RGB")
    if is_dark(im0):
        from PIL import ImageOps
        image = os.path.join(tmp, "inverted.png")
        ImageOps.invert(im0).save(image)
    passes = [("default", image, base, 1)]
    if not quick:
        passes += [("adaptive", image, base + ["-i"], 1), ("jaggy", image, base + ["-j"], 1)]
        from PIL import ImageFilter

        im = Image.open(image).convert("RGB")
        dn = os.path.join(tmp, "denoise.png")
        g = im.convert("L").filter(ImageFilter.MedianFilter(3))
        g.point(lambda v: 0 if v < 170 else 255).save(dn)
        passes.append(("denoise", dn, base, 1))
        if max(im.size) < 2500:
            up = os.path.join(tmp, "up2x.png")
            im.resize((im.width * 2, im.height * 2), Image.LANCZOS).save(up)
            passes.append(("upscale2x", up, base, 2))
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=min(len(passes), os.cpu_count() or 2)) as ex:
        runs = list(ex.map(lambda ps: (ps, run_osra(ps[1], ps[2])), passes))
    # A group is one reading of one drawing: same InChIKey AND same place.
    # (Two different drawings of the same compound stay separate.)
    hits, invalid = [], []
    for (label, path, extra, scale), found in runs:
        for h in found:
            if "box_px" in h and scale != 1:
                h["box_px"] = [v // scale for v in h["box_px"]]
            info = summarize(h["smiles"], brief=True)
            if not info["valid"]:
                invalid.append({**h, **info, "pass": label})
                continue
            box = h.get("box_px")
            g = next((g for g in hits if g["inchikey"] == info["inchikey"]
                      and ((box is None and g.get("box_px") is None)
                           or (box and g.get("box_px") and _overlap(box, g["box_px"]) > 0.5))), None)
            if g is None:
                g = {**h, **info, "passes": []}
                hits.append(g)
            if label not in g["passes"]:
                g["passes"].append(label)
            if h.get("confidence", -1e9) > g.get("confidence", -1e9):
                g.update({k: h[k] for k in ("confidence", "box_px", "bond_length_px") if k in h})
    for g in hits:
        g["votes"] = len(g["passes"])
        if "box_px" not in g:
            g["warning"] = "no location reported; cannot be tied to a label"
    for g in hits:
        if "*" in g["canonical_smiles"]:
            g["warning"] = "contains * (unread atom: label, compound number, or cropped edge)"
        else:
            # all-carbon readings with no aromatic ring are a classic OSRA failure
            mol = Chem.MolFromSmiles(g["canonical_smiles"])
            if mol and all(at.GetSymbol() == "C" for at in mol.GetAtoms()) \
                    and rdMolDescriptors.CalcNumAromaticRings(mol) == 0 and mol.GetNumHeavyAtoms() > 6:
                g["warning"] = "hydrocarbon with no aromatic ring: often a misread of an arene or a text label"
    hits.sort(key=lambda g: ("*" in g["canonical_smiles"], -g["votes"], -g.get("confidence", 0)))
    hits = merge_by_location(hits)
    if not hits:
        hits = invalid[:3]
    return hits


def _overlap(a, b):
    """Intersection over the smaller box (1.0 when one box contains the other)."""
    ix = max(0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0, min(a[3], b[3]) - max(a[1], b[1]))
    small = min((a[2] - a[0]) * (a[3] - a[1]), (b[2] - b[0]) * (b[3] - b[1]))
    return ix * iy / small if small > 0 else 0


def _iou(a, b):
    ix = max(0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0, min(a[3], b[3]) - max(a[1], b[1]))
    inter = ix * iy
    union = (a[2] - a[0]) * (a[3] - a[1]) + (b[2] - b[0]) * (b[3] - b[1]) - inter
    return inter / union if union else 0


def merge_by_location(hits, thresh=0.4):
    """Different passes can read the same drawing differently; keep one
    entry per drawing (the best-ranked reading) and list the rest as
    alternatives. Output is in reading order (top-to-bottom, left-to-right)."""
    merged = []
    for h in hits:  # already best-first
        box = h.get("box_px")
        home = next((m for m in merged if box and m.get("box_px")
                     and (_iou(box, m["box_px"]) > thresh or _overlap(box, m["box_px"]) > 0.85)), None)
        if home is None:
            merged.append(h)
        else:
            home.setdefault("alternatives", []).append(
                {k: h[k] for k in ("canonical_smiles", "votes", "confidence", "passes") if k in h})
    if merged and all("box_px" in m for m in merged):
        # Cluster into rows by vertical centre, then read each row left to right.
        rows = []
        for m in sorted(merged, key=lambda m: (m["box_px"][1] + m["box_px"][3]) / 2):
            cy, h = (m["box_px"][1] + m["box_px"][3]) / 2, m["box_px"][3] - m["box_px"][1]
            if rows and abs(cy - rows[-1][0]) < 0.5 * max(h, rows[-1][1]):
                rows[-1][2].append(m)
            else:
                rows.append([cy, h, [m]])
        merged = [m for _, _, row in rows for m in sorted(row, key=lambda m: m["box_px"][0])]
    return merged


def run_molscribe(image):
    import torch
    from molscribe import MolScribe

    ckpt = os.environ.get("MOLSCRIBE_CKPT")
    if not ckpt:
        from huggingface_hub import hf_hub_download
        ckpt = hf_hub_download("yujieq/MolScribe", "swin_base_char_aux_1m680k.pth")
    model = MolScribe(ckpt, device=torch.device("cuda" if torch.cuda.is_available() else "cpu"))
    r = model.predict_image_file(image, return_confidence=True)
    return [{"engine": "molscribe", "smiles": r.get("smiles"),
             "confidence": r.get("confidence"), "molfile": r.get("molfile")}]


def run_decimer(image):
    from DECIMER import predict_SMILES

    return [{"engine": "decimer", "smiles": predict_SMILES(image)}]


def cmd_recognize(a):
    avail = engines_available()
    wanted = a.engine.split(",") if a.engine != "all" else [e for e in avail if e != "rdkit"]
    results, notes = [], []
    for e in wanted:
        if not avail.get(e):
            notes.append(f"{e}: not installed")
            continue
        try:
            if e == "osra":
                hits = osra_ensemble(a.image, a.adaptive, a.jaggy, a.rotate, a.quick)
            elif e == "molscribe":
                hits = run_molscribe(a.image)
            else:
                hits = run_decimer(a.image)
            for h in hits:
                if "valid" not in h:
                    h.update(summarize(h.get("smiles") or "", brief=True))
            results.extend(hits)
            if not hits:
                notes.append(f"{e}: no structure found")
        except Exception as ex:
            notes.append(f"{e}: failed ({type(ex).__name__}: {ex})")
    from PIL import Image

    with Image.open(a.image) as im:
        size = list(im.size)
    out({"image": a.image, "image_size": size, "engines": avail, "results": results, "notes": notes})


# ---------------------------------------------------------------------- check

GROUPS = {
    "carboxylic acid": "[CX3](=O)[OX2H1]",
    "ester": "[#6][CX3](=O)[OX2H0][#6]",
    "lactone": "[#6;R][CX3;R](=O)[OX2;R][#6;R]",
    "amide": "[NX3][CX3](=[OX1])[#6]",
    "lactam": "[NX3;R][CX3;R](=[OX1])",
    "carbamate": "[NX3][CX3](=O)[OX2][#6]",
    "urea": "[NX3][CX3](=O)[NX3]",
    "aldehyde": "[CX3H1](=O)[#6]",
    "ketone": "[#6][CX3](=O)[#6]",
    "alcohol": "[OX2H][CX4]",
    "phenol": "[OX2H]c",
    "ether": "[OD2]([#6])[#6]",
    "primary amine": "[NX3;H2][CX4]",
    "secondary amine": "[NX3;H1]([CX4])[CX4]",
    "tertiary amine": "[NX3;H0;!$(NC=O);!$(N-a)]([CX4])([CX4])[CX4]",
    "aniline": "[NX3;!$(NC=O)]c",
    "nitrile": "C#N",
    "nitro": "[N+](=O)[O-]",
    "azide": "N=[N+]=[N-]",
    "alkene": "[CX3]=[CX3]",
    "alkyne": "C#C",
    "aromatic ring": "a1aaaaa1",
    "heteroaromatic ring": "[a;!c]",
    "halide (F)": "[F]",
    "halide (Cl)": "[Cl]",
    "halide (Br)": "[Br]",
    "halide (I)": "[I]",
    "sulfonamide": "S(=O)(=O)N",
    "sulfone": "[#6]S(=O)(=O)[#6]",
    "thiol": "[SX2H]",
    "thioether": "[SX2]([#6])[#6]",
    "epoxide": "C1OC1",
    "acetal": "[CX4]([OX2][#6])[OX2][#6]",
    "boronic acid/ester": "B(O)O",
    "phosphate": "P(=O)(O)O",
    "silyl ether": "[Si]O",
    "imine": "[CX3]=[NX2]",
    "acyl halide": "C(=O)[Cl,Br,F,I]",
    "anhydride": "C(=O)OC(=O)",
}
_GROUP_PATTERNS = {k: Chem.MolFromSmarts(v) for k, v in GROUPS.items()} if HAVE_RDKIT else {}


def summarize(smiles, brief=False):
    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles, sanitize=False) if smiles else None
    if mol is None:
        return {"valid": False, "error": "could not parse SMILES"}
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return {"valid": False, "error": f"chemistry problem: {e}"}
    frags = Chem.GetMolFrags(mol, asMols=True)
    info = {
        "valid": True,
        "canonical_smiles": Chem.MolToSmiles(mol),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": round(Descriptors.MolWt(mol), 2),
        "inchikey": Chem.MolToInchiKey(mol) or None,
    }
    if brief:
        return info
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    info.update({
        "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
        "inchi": Chem.MolToInchi(mol) or None,
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "fragments": len(frags),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "formal_charge": Chem.GetFormalCharge(mol),
        "stereocenters": [{"atom": i + 1, "element": mol.GetAtomWithIdx(i).GetSymbol(), "label": l}
                          for i, l in centers],  # 1-based, like build
        "unassigned_stereocenters": sum(1 for _, l in centers if l == "?"),
        "functional_groups": sorted(k for k, p in _GROUP_PATTERNS.items()
                                    if mol.HasSubstructMatch(p)),
        "warnings": sanity_warnings(mol, frags),
    })
    return info


def sanity_warnings(mol, frags):
    """Cheap heuristics that often flag a misread structure."""
    w = []
    if len(frags) > 1:
        w.append(f"{len(frags)} disconnected fragments (salt/counter-ion, or a missed bond?)")
    for at in mol.GetAtoms():
        if at.GetNumRadicalElectrons():
            w.append(f"radical on atom {at.GetIdx()} ({at.GetSymbol()}): likely a missing H or bond")
        if at.GetSymbol() not in ("C", "H", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se", "Na", "K", "Li", "Mg", "Zn", "Cu", "Pd", "Sn", "Fe"):
            w.append(f"unusual element {at.GetSymbol()} at atom {at.GetIdx()}: check OCR of an atom label")
        if at.GetFormalCharge() and at.GetSymbol() == "C":
            w.append(f"charged carbon at atom {at.GetIdx()}")
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) > 8 and not any(len(r) <= 8 for r in ri.AtomRings() if set(r) & set(ring)):
            w.append(f"large isolated {len(ring)}-membered ring: check whether it is real or a misread")
    return sorted(set(w))


def cmd_check(a):
    smiles = a.smiles or [l for l in sys.stdin.read().split() if l]
    res = [{"input": s, **summarize(s)} for s in smiles]
    out(res if len(res) > 1 else res[0])


# ------------------------------------------------------------- render/compare

def draw(mol, path, size, legend=""):
    has_layout = mol.GetNumConformers() > 0
    if not has_layout:
        AllChem.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    o = d.drawOptions()
    if has_layout:
        o.useMolBlockWedging = True  # keep the wedges where the drawing has them
    o.addStereoAnnotation = True
    o.legendFontSize = 18
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, legend=legend)
    d.FinishDrawing()
    with open(path, "wb") as f:
        f.write(d.GetDrawingText())


def cmd_render(a):
    from PIL import Image, ImageDraw

    mols = []
    for i, s in enumerate(a.smiles):
        m = Chem.MolFromMolFile(s) if s.endswith(".mol") and os.path.exists(s) else Chem.MolFromSmiles(s)
        if m is None:
            out({"error": f"invalid SMILES or molfile: {s}"})
            sys.exit(1)
        mols.append(m)
    w, h = a.size
    tmpdir = tempfile.mkdtemp()
    tiles = []
    for i, m in enumerate(mols):
        p = os.path.join(tmpdir, f"{i}.png")
        legend = a.labels[i] if a.labels and i < len(a.labels) else rdMolDescriptors.CalcMolFormula(m)
        draw(m, p, (w, h), legend)
        tiles.append(Image.open(p).convert("RGB"))
    if a.original:
        src = Image.open(a.original).convert("RGB")
        if is_dark(src):
            from PIL import ImageOps
            src = ImageOps.invert(src)
        src.thumbnail((w, h))  # keeps aspect ratio
        pad = Image.new("RGB", (w, h), "white")
        pad.paste(src, ((w - src.width) // 2, (h - src.height) // 2))
        ImageDraw.Draw(pad).text((8, 8), "original", fill=(200, 0, 0))
        tiles.insert(0, pad)
    cols = min(len(tiles), a.cols)
    rows = (len(tiles) + cols - 1) // cols
    canvas = Image.new("RGB", (w * cols, h * rows), "white")
    dr = ImageDraw.Draw(canvas)
    for i, t in enumerate(tiles):
        x, y = (i % cols) * w, (i // cols) * h
        canvas.paste(t, (x, y))
        dr.rectangle([x, y, x + w - 1, y + h - 1], outline=(190, 190, 190), width=1)
    canvas.save(a.out)
    out({"image": a.out, "panels": (["original"] if a.original else []) + list(a.smiles)})


_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048) if HAVE_RDKIT else None


def fp(m):
    return _MORGAN.GetFingerprint(m)


def is_mirror(mol, other_key):
    """True if inverting every stereocentre of mol gives the other InChIKey
    (same compound, opposite absolute configuration: a match for racemic or
    relative-stereo drawings)."""
    m = Chem.RWMol(mol)
    for at in m.GetAtoms():
        t = at.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return Chem.MolToInchiKey(m) == other_key


def cmd_compare(a):
    items = []
    for s in a.smiles:
        m = Chem.MolFromSmiles(s)
        items.append({"smiles": s, "mol": m,
                      "inchikey": Chem.MolToInchiKey(m) if m else None,
                      "formula": rdMolDescriptors.CalcMolFormula(m) if m else None})
    pairs = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            A, B = items[i], items[j]
            if not (A["mol"] and B["mol"]):
                pairs.append({"a": i, "b": j, "error": "invalid SMILES"})
                continue
            ka, kb = A["inchikey"], B["inchikey"]
            pairs.append({
                "a": i, "b": j,
                "identical": ka == kb,
                "same_skeleton_ignoring_stereo": ka.split("-")[0] == kb.split("-")[0],
                "same_formula": A["formula"] == B["formula"],
                "tanimoto": round(DataStructs.TanimotoSimilarity(fp(A["mol"]), fp(B["mol"])), 3),
                "enantiomers": ka != kb and is_mirror(A["mol"], kb),
            })
    out({"candidates": [{k: v for k, v in it.items() if k != "mol"} for it in items],
         "pairs": pairs})


# ---------------------------------------------------------------------- build

def cmd_build(a):
    """Build a molecule from atoms placed at their image coordinates plus
    bonds with wedge/hash marks, so RDKit assigns R/S from the drawing itself
    instead of from a hand-written SMILES. Also keeps the drawing's layout,
    so `render` of the saved molfile looks like the original."""
    spec = json.load(open(a.spec) if a.spec != "-" else sys.stdin)
    atoms, bonds = spec["atoms"], spec["bonds"]
    stereo_code = {None: 0, "": 0, "wedge": 1, "hash": 6, "wavy": 4}
    lines = ["", "  chemvision", "",
             f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000"]
    rlabels = {}
    for n, at in enumerate(atoms, 1):
        sym, x, y = at[0], float(at[1]), float(at[2])
        if re.fullmatch(r"R\d*|X|Ar|\*", sym):  # R-group / attachment point
            rlabels[n] = sym
            sym = "R#" if sym != "*" else "*"
        charge = at[3] if len(at) > 3 else 0
        lines.append(f"{x / 50:10.4f}{-y / 50:10.4f}{0:10.4f} {sym:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
        if charge:
            lines.append(None)  # placeholder, charges go in M  CHG below
    lines = [l for l in lines if l is not None]
    for b in bonds:
        i, j, order = int(b[0]), int(b[1]), int(b[2])
        st = stereo_code[b[3] if len(b) > 3 else None]
        lines.append(f"{i:3d}{j:3d}{order:3d}{st:3d}")
    for idx, at in enumerate(atoms, 1):
        if len(at) > 3 and at[3]:
            lines.append(f"M  CHG  1 {idx:3d} {int(at[3]):3d}")
    lines.append("M  END")
    for n, lab in rlabels.items():
        num = re.sub(r"\D", "", lab)
        if num:
            lines.insert(-1, f"M  RGP  1 {n:3d} {int(num):3d}")
    m = Chem.MolFromMolBlock("\n".join(lines), removeHs=False)
    if m is not None:
        for n, lab in rlabels.items():  # keep R numbering visible in SMILES as [1*], [2*]
            at = m.GetAtomWithIdx(n - 1)
            num = re.sub(r"\D", "", lab)
            at.SetIsotope(int(num) if num else 0)
            at.SetProp("dummyLabel", lab)
    if m is None:
        out({"valid": False, "error": "RDKit rejected the atom/bond spec (check valences and indices)"})
        sys.exit(1)
    m = Chem.RemoveHs(m)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    if a.out:
        Chem.MolToMolFile(m, a.out)
    smi = Chem.MolToSmiles(m)
    res = {"smiles": smi, **summarize(smi)}
    if a.out:
        res["molfile"] = a.out
    out(res)


# ----------------------------------------------------------- annotate/export

def cmd_annotate(a):
    import annotations as ann

    if a.pair:
        recs = []
        for struct, text in a.pair:
            r = ann.annotate_text(text)
            ids = r.pop("ids")
            rec = {"id": ids[0] if ids else None, **{k: v for k, v in r.items()}}
            if struct.endswith(".mol") and os.path.exists(struct) and not HAVE_RDKIT:
                rec.update({"smiles": None, "molfile": struct, "note": "RDKit missing: SMILES not derived"})
            elif struct.endswith(".mol") and os.path.exists(struct):
                m = Chem.MolFromMolFile(struct)
                rec.update({"smiles": Chem.MolToSmiles(m) if m else None, "molfile": struct})
            else:
                rec["smiles"] = struct
            rec["source_text"] = rec.pop("text")
            recs.append(rec)
        out(recs)
        return
    if a.text:
        res = [ann.annotate_text(t) for t in a.text]
        out(res if len(res) > 1 else res[0])
        return
    if not a.pdf:
        out({"error": "give --text, --pair, or --pdf"})
        sys.exit(1)
    import pymupdf as fitz

    doc = fitz.open(a.pdf)
    page = doc[a.page - 1]
    W, H = page.rect.width, page.rect.height
    has_text = bool(page.get_text().strip())
    if not has_text and not a.labels:
        out({"note": ("No text layer on this page (scanned, or the figure is a raster image). "
                      "Read the labels visually, then either pass them with their positions via "
                      "--labels, or pair them with your structures via --pair."),
             "pdf": a.pdf, "page": a.page})
        return

    boxes, img_size, rec = [], None, None
    if a.from_recognize:
        rec = json.load(open(a.from_recognize))
        boxes = [h["box_px"] for h in rec.get("results", []) if "box_px" in h]
        img_size = rec.get("image_size")
        units = "px"
    elif a.boxes:
        boxes = json.loads(a.boxes)
        units = a.units

    def to_pt(b):
        if units == "pt":
            return list(b)
        if units == "frac":
            return [b[0] * W, b[1] * H, b[2] * W, b[3] * H]
        if a.image_rect:  # px in an embedded image placed at image_rect on the page
            x0, y0, x1, y1 = [float(v) for v in a.image_rect.split(",")]
            iw, ih = (img_size or a.image_size)
            sx, sy = (x1 - x0) / iw, (y1 - y0) / ih
            return [x0 + b[0] * sx, y0 + b[1] * sy, x0 + b[2] * sx, y0 + b[3] * sy]
        f = 72.0 / a.dpi  # px in a full-page render at --dpi
        return [v * f for v in b]

    if a.image_rect and not (img_size or a.image_size):
        out({"error": "--image-rect needs the image size: use --from-recognize, or pass --image-size W H"})
        sys.exit(1)
    pt_boxes = [to_pt(b) for b in boxes]
    if a.labels:
        labels = json.loads(a.labels)
        for lab in labels:
            lab["box"] = to_pt(lab["box"])
        res = ann.annotate_labels(pt_boxes, labels, a.max_dist) if pt_boxes else {"lines": labels}
    else:
        res = ann.annotate_page(page, pt_boxes or None, a.max_dist)
    for s in res.get("structures", []):
        s["box_pt"] = [round(v, 1) for v in s.pop("box")]
    res = {"pdf": a.pdf, "page": a.page, "page_size_pt": [W, H], **res}
    if rec:
        for s, h in zip(res.get("structures", []), [h for h in rec["results"] if "box_px" in h]):
            s.update({k: h[k] for k in ("smiles", "canonical_smiles", "formula", "inchikey", "votes") if k in h})
    out(res)


EXPORT_CORE = ["id", "smiles", "name", "formula", "mw", "inchikey",
               "yield", "yield_assumed", "yield_nmr", "conversion", "ee", "er", "ee_from_er", "dr", "de", "rr", "ez",
               "conditions", "temperature", "time", "scale", "loading",
               "page", "figure", "confidence", "source_text", "notes"]


def cmd_export(a):
    """Write compound records (JSON list) to CSV, SDF or JSON, filling in
    formula/MW/InChIKey from the SMILES and flagging invalid ones."""
    recs = json.load(open(a.records) if a.records != "-" else sys.stdin)
    if isinstance(recs, dict):
        recs = recs.get("compounds") or recs.get("structures") or [recs]
    rows, problems = [], []
    for r in recs:
        r = dict(r)
        r.setdefault("id", r.pop("compound_id", None))
        smi = r.get("smiles") or r.get("canonical_smiles")
        if smi:
            info = summarize(smi, brief=True)
            if info["valid"]:
                r.update({"smiles": info["canonical_smiles"], "formula": info["formula"],
                          "mw": info["mw"], "inchikey": info["inchikey"]})
                if "*" in info["canonical_smiles"]:
                    r.update({"markush": True, "mw": None, "formula": info["formula"] + " (core only)"})
                m = Chem.MolFromSmiles(info["canonical_smiles"])
                n_stereo = len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False))
                if r.get("dr") and n_stereo < 2:
                    problems.append(f"{r.get('id')}: dr {r['dr']} given but the structure has {n_stereo} stereocentre(s)")
                if (r.get("ee") or r.get("er")) and n_stereo == 0:
                    problems.append(f"{r.get('id')}: ee/er given but the structure has no stereocentres")
            else:
                problems.append(f"{r.get('id')}: {info['error']}")
        r.pop("canonical_smiles", None)
        rows.append(r)
    cols = [c for c in EXPORT_CORE if any(c in r for r in rows)]
    cols += sorted({k for r in rows for k in r} - set(cols) - {"molfile", "bbox", "box_pt", "box_px", "votes", "passes", "alternatives", "valid"})
    fmt = a.format or os.path.splitext(a.out)[1].lstrip(".").lower()
    if fmt == "csv":
        import csv
        with open(a.out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow({k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.items()})
    elif fmt == "sdf":
        w = Chem.SDWriter(a.out)
        for r in rows:
            mf = r.get("molfile")
            m = Chem.MolFromMolFile(mf) if mf and os.path.exists(mf) else None
            if m is None:  # no drawing layout available: lay out from SMILES
                m = Chem.MolFromSmiles(r.get("smiles") or "")
                if m is None:
                    continue
                AllChem.Compute2DCoords(m)
            if r.get("id"):
                m.SetProp("_Name", str(r["id"]))
            for k in cols:
                if k not in ("smiles",) and r.get(k) not in (None, ""):
                    m.SetProp(k, json.dumps(r[k]) if isinstance(r[k], (list, dict)) else str(r[k]))
            w.write(m)
        w.close()
    else:
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(rows, f, indent=2, ensure_ascii=False)
    out({"written": a.out, "format": fmt, "compounds": len(rows), "columns": cols, "problems": problems})


# ------------------------------------------------------------------- reaction

def cmd_reaction(a):
    from rdkit.Chem import rdChemReactions

    text = a.rxn.strip()
    try:
        rxn = rdChemReactions.ReactionFromSmarts(text, useSmiles=True)
    except Exception as e:
        out({"valid": False, "error": f"could not parse reaction: {e}"})
        sys.exit(1)
    sides = {"reactants": rxn.GetReactants(), "agents": rxn.GetAgents(), "products": rxn.GetProducts()}
    parts, bad = {}, []
    for side, mols in sides.items():
        parts[side] = []
        for m in mols:
            s = Chem.MolToSmiles(m)
            info = summarize(s, brief=True)
            if not info["valid"]:
                bad.append(f"{side}: {s}")
            parts[side].append({"smiles": s, **info})

    def atoms(mols):
        cnt = {}
        for m in mols:
            mm = Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(m)))
            for at in mm.GetAtoms():
                cnt[at.GetSymbol()] = cnt.get(at.GetSymbol(), 0) + 1
        return cnt

    r_atoms, p_atoms = atoms(sides["reactants"]), atoms(sides["products"][:1])
    diff = {el: p_atoms.get(el, 0) - r_atoms.get(el, 0)
            for el in sorted(set(r_atoms) | set(p_atoms))
            if p_atoms.get(el, 0) != r_atoms.get(el, 0)}
    result = {
        **{k: v for k, v in (("label", a.label), ("conditions", a.conditions), ("yield", a.yld)) if v},
        "valid": not bad,
        "invalid_components": bad,
        **parts,
        "balanced": not diff,
        "product_minus_reactant_atoms": diff,
        "diff_basis": ("first product only (list the main product first; others are treated as "
                       "co-products or isomers)" if len(sides["products"]) > 1 else "single product"),
        "note": ("Unbalanced is normal for schemes that omit by-products or reagents; "
                 "use the atom difference to sanity-check what the arrow implies."),
    }
    if a.out:
        img = Draw.ReactionToImage(rxn, subImgSize=(300, 250))
        img.save(a.out)
        result["image"] = a.out
    out(result)


# --------------------------------------------------------------------- lookup

def cmd_lookup(a):
    import urllib.parse
    import urllib.request

    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    q = a.query.strip()
    if re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", q):
        url = f"{base}/inchikey/{q}"
    elif a.smiles:
        url = f"{base}/smiles/{urllib.parse.quote(q, safe='')}"
    else:
        url = f"{base}/name/{urllib.parse.quote(q, safe='')}"
    url += "/property/IUPACName,Title,MolecularFormula,CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            data = json.load(r)
        out({"query": q, "pubchem": data["PropertyTable"]["Properties"]})
    except Exception as e:
        out({"query": q, "error": f"PubChem lookup failed: {e}"})


def cmd_engines(a):
    eng = engines_available()
    try:
        import pymupdf  # noqa: F401
        eng["pymupdf"] = True
    except Exception:
        eng["pymupdf"] = False
    eng["poppler"] = bool(shutil.which("pdftoppm"))
    eng["libreoffice"] = bool(shutil.which("soffice") or shutil.which("libreoffice"))
    if not HAVE_RDKIT:
        eng["rdkit_missing"] = {**RDKIT_MISSING, "detail": RDKIT_ERROR}
    if not eng["pymupdf"]:
        eng["pymupdf_note"] = ("PyMuPDF missing: `pages` falls back to Poppler (pdftoppm/pdftotext/pdfimages); "
                               "`annotate --pdf` needs PyMuPDF.") if eng["poppler"] else \
            "PyMuPDF and Poppler both missing: PDF pages cannot be rendered."
    out(eng)


# ------------------------------------------------- pages fallback / grid overlay

def pages_poppler(a):
    """`pages` without PyMuPDF: render with Poppler's command-line tools."""
    if not shutil.which("pdftoppm"):
        out({"error": "Neither PyMuPDF nor Poppler (pdftoppm) is available; cannot render PDF pages."})
        sys.exit(3)
    info = subprocess.run(["pdfinfo", a.pdf], capture_output=True, text=True).stdout
    m = re.search(r"^Pages:\s+(\d+)", info, re.M)
    n = int(m.group(1)) if m else 1
    os.makedirs(a.outdir, exist_ok=True)
    pages = parse_range(a.pages, n) if a.pages else range(n)
    stem = os.path.splitext(os.path.basename(a.pdf))[0]
    result = []
    from PIL import Image
    for i in pages:
        base = os.path.join(a.outdir, f"{stem}_p{i + 1}")
        subprocess.run(["pdftoppm", "-png", "-r", str(a.dpi), "-f", str(i + 1), "-l", str(i + 1),
                        "-singlefile", a.pdf, base], check=True)
        path = base + ".png"
        with Image.open(path) as im:
            entry = {"page": i + 1, "png": path, "size": list(im.size)}
        if a.grid:
            entry["grid"] = draw_grid(path, base + "_grid.png")
        if a.images and shutil.which("pdfimages"):
            idir = base + "_img"
            os.makedirs(idir, exist_ok=True)
            subprocess.run(["pdfimages", "-png", "-f", str(i + 1), "-l", str(i + 1), a.pdf,
                            os.path.join(idir, "img")], check=False)
            imgs = []
            for f in sorted(os.listdir(idir)):
                fp_ = os.path.join(idir, f)
                with Image.open(fp_) as im:
                    if im.width >= 80 and im.height >= 80:
                        imgs.append({"png": fp_, "size": list(im.size), "page_rect_pt": None})
            entry["embedded_images"] = imgs
        if a.text and shutil.which("pdftotext"):
            entry["text"] = subprocess.run(["pdftotext", "-layout", "-f", str(i + 1), "-l", str(i + 1),
                                            a.pdf, "-"], capture_output=True, text=True).stdout
        result.append(entry)
    out({"pdf": a.pdf, "page_count": n, "backend": "poppler (PyMuPDF not installed)", "pages": result,
         "note": "embedded image page positions (page_rect_pt) are unavailable without PyMuPDF"})


def draw_grid(src, dst, step=0.1, max_w=1400):
    """Write a downscaled copy of `src` with labeled gridlines every `step`
    (as a fraction of width/height). Labels give the fraction and the
    full-resolution pixel coordinate, so a box chosen on this view can be passed
    to `crop --box` directly as fractions (preferred) or full-res pixels."""
    from PIL import Image, ImageDraw, ImageFont

    im = Image.open(src).convert("RGB")
    W, H = im.size
    k = min(1.0, max_w / W)
    v = im.resize((int(W * k), int(H * k)), Image.LANCZOS) if k < 1 else im.copy()
    w, h = v.size
    d = ImageDraw.Draw(v)
    try:
        font = ImageFont.truetype("DejaVuSans.ttf", max(10, w // 90))
    except Exception:
        font = ImageFont.load_default()
    n = int(round(1 / step))
    for i in range(n + 1):
        f = i * step
        x, y = int(f * (w - 1)), int(f * (h - 1))
        col = (220, 0, 0) if i % 5 == 0 else (255, 120, 120)
        d.line([(x, 0), (x, h)], fill=col, width=1)
        d.line([(0, y), (w, y)], fill=col, width=1)
        if 0 < i < n:
            d.text((x + 2, 2), f"{f:.1f}|{int(f * W)}", fill=(200, 0, 0), font=font)
            d.text((2, y + 2), f"{f:.1f}|{int(f * H)}", fill=(200, 0, 0), font=font)
    v.save(dst)
    return {"grid_png": dst, "full_size": [W, H], "view_size": [w, h], "step": step,
            "how_to_use": "labels read fraction|full-res-px; pass crop --box as fractions, e.g. 0.30,0.40,0.55,0.62"}


def cmd_grid(a):
    dst = a.out or os.path.splitext(a.image)[0] + "_grid.png"
    out(draw_grid(a.image, dst, step=a.step, max_w=a.max_width))


# ---------------------------------------------------------------- docx figures

def cmd_docx_figures(a):
    """Extract every embedded figure from a .docx in document order and convert
    vector formats (EMF/WMF, typical for pasted ChemDraw) to PNG at --dpi.
    Each EMF is converted via PDF so any real text it carries can be reported."""
    import zipfile
    import xml.etree.ElementTree as ET
    from PIL import Image

    os.makedirs(a.outdir, exist_ok=True)
    z = zipfile.ZipFile(a.docx)
    rels = {}
    if "word/_rels/document.xml.rels" in z.namelist():
        for r in ET.fromstring(z.read("word/_rels/document.xml.rels")):
            rels[r.get("Id")] = r.get("Target")
    W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
    order, context = [], {}
    body = ET.fromstring(z.read("word/document.xml"))
    paras = list(body.iter(W + "p"))
    for pi, para in enumerate(paras):
        for el in para.iter():
            rid = el.get(R + "embed") or el.get(R + "id")
            if rid and rid in rels and "media/" in rels[rid] and rid not in order:
                order.append(rid)
                txt = lambda q: "".join(t.text or "" for t in q.iter(W + "t")).strip()
                near = [txt(paras[j]) for j in (pi - 1, pi, pi + 1) if 0 <= j < len(paras)]
                context[rid] = " … ".join(x for x in near if x)[:300]
    soffice = shutil.which("soffice") or shutil.which("libreoffice")
    figs = []
    for n, rid in enumerate(order, 1):
        target = rels[rid]
        name = "word/" + target if not target.startswith("word/") else target
        ext = os.path.splitext(target)[1].lower()
        raw = os.path.join(a.outdir, f"fig{n}{ext}")
        with open(raw, "wb") as fh:
            fh.write(z.read(name))
        rec = {"figure": n, "source": target, "format": ext.lstrip("."), "nearby_text": context.get(rid, "")}
        try:
            if ext in (".png", ".jpg", ".jpeg", ".gif", ".bmp", ".tif", ".tiff"):
                png = os.path.join(a.outdir, f"fig{n}.png")
                Image.open(raw).convert("RGB").save(png)
                rec["png"] = png
            elif ext in (".emf", ".wmf", ".svg") and soffice:
                tmp = tempfile.mkdtemp()
                subprocess.run([soffice, "--headless", "--convert-to", "pdf", "--outdir", tmp, raw],
                               capture_output=True, timeout=300)
                base = os.path.join(a.outdir, f"fig{n}")
                pdf = base + ".pdf"
                shutil.move(os.path.join(tmp, f"fig{n}.pdf"), pdf)
                subprocess.run(["pdftoppm", "-png", "-r", str(a.dpi), "-singlefile", pdf, base], check=True)
                rec["png"] = base + ".png"
                rec["pdf"] = pdf
                if shutil.which("pdftotext"):
                    t = subprocess.run(["pdftotext", pdf, "-"], capture_output=True, text=True).stdout.strip()
                    rec["text_layer"] = bool(t)
                    if t:
                        rec["title"] = t.splitlines()[0].strip()[:120]
                        rec["text"] = t[:6000]
                        rec["note"] = ("real text: labels/values are exact. `annotate --pdf <this pdf> --page 1` "
                                       "ties them to structures (needs PyMuPDF); otherwise search `text`")
                    else:
                        rec["note"] = "no real text in this figure: read labels and values visually"
            else:
                rec["error"] = f"cannot convert {ext} (LibreOffice missing?)"
            if "png" in rec:
                from PIL import ImageChops
                im = Image.open(rec["png"]).convert("RGB")
                bg = Image.new("RGB", im.size, (255, 255, 255))
                bb = ImageChops.difference(im, bg).getbbox()
                if bb and (bb[2] - bb[0]) * (bb[3] - bb[1]) < 0.95 * im.width * im.height:
                    pad = 12
                    bb = (max(0, bb[0] - pad), max(0, bb[1] - pad), min(im.width, bb[2] + pad), min(im.height, bb[3] + pad))
                    im.crop(bb).save(rec["png"])
                    rec["trimmed_from"] = list(im.size)
                    im = Image.open(rec["png"])
                rec["size"] = list(im.size)
                if a.grid:
                    rec["grid"] = draw_grid(rec["png"], rec["png"][:-4] + "_grid.png")
        except Exception as e:
            rec["error"] = f"{type(e).__name__}: {e}"
        figs.append(rec)
    out({"docx": a.docx, "figure_count": len(figs), "figures": figs,
         "hint": "If the document's ChemDraw source (.cdxml/.cdx) is available, run `cdxml` on it instead: exact text and structures."})


# ----------------------------------------------------------------------- cdxml

def cmd_cdxml(a):
    """Read a ChemDraw CDXML file: every structure (fragment) with its bounding
    box, the caption text objects (compound labels, yields, ee ...), and the
    labels tied to structures by position. With RDKit, also SMILES per fragment."""
    import xml.etree.ElementTree as ET
    import annotations as ann

    root = ET.parse(a.cdxml).getroot()
    parent = {c: p for p in root.iter() for c in p}

    def bbox(el):
        b = el.get("BoundingBox")
        return [float(v) for v in b.split()] if b else None

    def inside_node(el):
        q = parent.get(el)
        while q is not None:
            if q.tag == "n":
                return True
            q = parent.get(q)
        return False

    frags = []
    for f in root.iter("fragment"):
        if inside_node(f):  # abbreviation expansion (Bn, Ts ...) inside an atom
            continue
        b = bbox(f)
        if b and f.find("n") is not None:
            frags.append({"fragment_id": f.get("id"), "box_pt": b,
                          "atoms": sum(1 for _ in f.findall("n"))})
    labels = []
    for t in root.iter("t"):
        if inside_node(t):  # atom labels (N, Me, OMe ...)
            continue
        text = "".join(s_.text or "" for s_ in t.iter("s")).strip()
        b = bbox(t)
        if text and b:
            labels.append({"text": text, "box": b})

    res = ann.annotate_labels([f["box_pt"] for f in frags], labels, a.max_dist) if frags else \
        {"structures": [], "unassigned_lines": labels}
    structs = []
    for f, s_ in zip(frags, res["structures"]):
        s_.pop("box", None)
        structs.append({**f, **s_})

    notes = []
    if HAVE_RDKIT:
        try:
            mols = Chem.MolsFromCDXMLFile(a.cdxml, sanitize=True, removeHs=True)
            by_id, unmapped = {}, []
            for m in mols:
                if m is None:
                    continue
                fid = m.GetProp("CDXML_FRAG_ID") if m.HasProp("CDXML_FRAG_ID") else None
                (by_id.__setitem__(str(fid), m) if fid is not None else unmapped.append(m))
            for s_ in structs:
                m = by_id.get(str(s_["fragment_id"]))
                if m is not None:
                    s_["smiles"] = Chem.MolToSmiles(m)
            if unmapped:
                notes.append(f"{len(unmapped)} RDKit molecules could not be tied to a fragment id; "
                             "their SMILES are listed under unmapped_smiles")
                res["unmapped_smiles"] = [Chem.MolToSmiles(m) for m in unmapped]
        except Exception as e:
            notes.append(f"RDKit could not parse structures: {type(e).__name__}: {e}")
    else:
        notes.append("RDKit missing: labels and values are exact, but no SMILES were derived. "
                     "Read the structures from the figure and use `annotate --pair`.")

    if a.id:
        want = set(a.id)
        structs = [s_ for s_ in structs if s_.get("compound_id") in want]
    labelled = [s_ for s_ in structs if s_.get("compound_id")]
    seen = {}
    for s_ in labelled:
        seen.setdefault(s_["compound_id"], []).append(s_["fragment_id"])
    dups = {k: v for k, v in seen.items() if len(v) > 1}
    if dups:
        notes.append("Some compound IDs label more than one structure (several schemes or old versions "
                     "on one canvas). Check box_pt positions to pick the right one: " +
                     ", ".join(f"{k} x{len(v)}" for k, v in sorted(dups.items())))
    out({"cdxml": a.cdxml, "fragments": len(frags), "labelled_structures": len(labelled),
         "structures": labelled if not a.all else structs,
         "unassigned_text": [] if a.id else [u.get("text") for u in res["unassigned_lines"]][:200],
         "notes": notes})


# ----------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)

    p = sp.add_parser("pages", help="render PDF pages to PNG")
    p.add_argument("pdf")
    p.add_argument("--outdir", default="chemvision_pages")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--pages", help="e.g. 1,3-5 (1-based)")
    p.add_argument("--images", action="store_true", help="also extract embedded raster images")
    p.add_argument("--text", action="store_true", help="include page text (for compound numbers, captions)")
    p.add_argument("--grid", action="store_true", help="also write a labeled grid overlay per page (for crop boxes)")
    p.set_defaults(func=cmd_pages)

    p = sp.add_parser("grid", help="draw a labeled coordinate grid on an image")
    p.add_argument("image")
    p.add_argument("--step", type=float, default=0.1, help="grid spacing as a fraction (default 0.1)")
    p.add_argument("--max-width", type=int, default=1400, help="width of the overlay view in px")
    p.add_argument("-o", "--out")
    p.set_defaults(func=cmd_grid)

    p = sp.add_parser("docx-figures", help="extract embedded figures from a .docx as PNG")
    p.add_argument("docx")
    p.add_argument("--outdir", default="chemvision_docx")
    p.add_argument("--dpi", type=int, default=300, help="render resolution for EMF/WMF figures")
    p.add_argument("--grid", action="store_true", help="also write a grid overlay per figure")
    p.set_defaults(func=cmd_docx_figures)

    p = sp.add_parser("cdxml", help="read labels (and SMILES with RDKit) from a ChemDraw .cdxml")
    p.add_argument("cdxml")
    p.add_argument("--id", nargs="+", help="only these compound IDs, e.g. --id 3m 4a")
    p.add_argument("--all", action="store_true", help="include unlabelled structures")
    p.add_argument("--max-dist", type=float, help="max label distance from a structure, in points")
    p.set_defaults(func=cmd_cdxml)

    p = sp.add_parser("crop", help="crop and upscale a region")
    p.add_argument("image")
    p.add_argument("--box", required=True, help="x0,y0,x1,y1 in pixels or 0-1 fractions")
    p.add_argument("--scale", type=float, default=2.0)
    p.add_argument("--pad", type=int, default=10, help="pixels added on every side (0 for an exact box)")
    p.add_argument("--invert", action="store_true", help="invert colours (white-on-black figures)")
    p.add_argument("--no-auto-invert", action="store_true", help="do not auto-invert dark crops")
    p.add_argument("-o", "--out", default="crop.png")
    p.set_defaults(func=cmd_crop)

    p = sp.add_parser("recognize", help="run OCSR engines on an image")
    p.add_argument("image")
    p.add_argument("--engine", default="all", help="osra, molscribe, decimer, comma list, or all")
    p.add_argument("--adaptive", action="store_true", help="OSRA adaptive thresholding (scans, low contrast)")
    p.add_argument("--jaggy", action="store_true", help="OSRA extra thinning (low quality documents)")
    p.add_argument("--rotate", type=float, default=0, help="rotate clockwise by N degrees before OSRA")
    p.add_argument("--quick", action="store_true", help="OSRA single pass instead of the preprocessing ensemble")
    p.set_defaults(func=cmd_recognize)

    p = sp.add_parser("check", help="validate and describe SMILES")
    p.add_argument("smiles", nargs="*")
    p.set_defaults(func=cmd_check)

    p = sp.add_parser("render", help="draw SMILES for visual comparison")
    p.add_argument("smiles", nargs="+", help="SMILES strings or .mol files (molfiles keep their layout)")
    p.add_argument("--original", help="source image to show as the first panel")
    p.add_argument("--labels", nargs="*")
    p.add_argument("--size", type=int, nargs=2, default=[450, 380])
    p.add_argument("--cols", type=int, default=3, help="panels per row")
    p.add_argument("-o", "--out", default="render.png")
    p.set_defaults(func=cmd_render)

    p = sp.add_parser("build", help="build a molecule from drawn atom positions + wedge/hash bonds")
    p.add_argument("spec", help='JSON file ("-" for stdin): {"atoms": [["C", x, y], ["N", x, y, +1], ...], '
                                '"bonds": [[1, 2, 1], [2, 3, 2], [3, 4, 1, "wedge"], ...]} (1-based, '
                                'wedge/hash narrow end at the first atom; explicit ["H", x, y] atoms allowed)')
    p.add_argument("-o", "--out", help="save a .mol file (keeps the drawing's layout for render)")
    p.set_defaults(func=cmd_build)

    p = sp.add_parser("compare", help="agreement between candidate SMILES")
    p.add_argument("smiles", nargs="+")
    p.set_defaults(func=cmd_compare)

    p = sp.add_parser("reaction", help="parse/check a reaction SMILES")
    p.add_argument("rxn", help="reactants>agents>products")
    p.add_argument("--conditions", help='free text, e.g. "Pd(PPh3)4, K2CO3, dioxane/H2O, 90 °C, 12 h"')
    p.add_argument("--yield", dest="yld", help="e.g. 87%%")
    p.add_argument("--label", help='step label, e.g. "1 + 2 -> 3"')
    p.add_argument("-o", "--out", help="write a scheme image")
    p.set_defaults(func=cmd_reaction)

    p = sp.add_parser("annotate", help="compound IDs, yields, ee/er, dr next to structures")
    p.add_argument("--text", nargs="+", help='label text read from a figure, e.g. "3a, 92%%, 95%% ee"')
    p.add_argument("--pair", nargs=2, action="append", metavar=("STRUCTURE", "TEXT"),
                   help='SMILES or .mol file plus its label text; repeatable; output is ready for export')
    p.add_argument("--labels", help='raster pages: JSON [{"text": "3a, 92%%", "box": [x0,y0,x1,y1]}, ...] '
                                    'read visually, boxes in the same units as the structure boxes')
    p.add_argument("--image-rect", help="x0,y0,x1,y1 (pt): where the image that px boxes came from sits on "
                                        "the page (page_rect_pt from `pages --images`)")
    p.add_argument("--image-size", type=int, nargs=2, help="W H of that image if not using --from-recognize")
    p.add_argument("--pdf", help="PDF with a text layer (vector figures)")
    p.add_argument("--page", type=int, default=1, help="1-based page")
    p.add_argument("--boxes", help="JSON list of structure boxes [[x0,y0,x1,y1], ...]")
    p.add_argument("--units", choices=["px", "pt", "frac"], default="px", help="units of --boxes")
    p.add_argument("--dpi", type=int, default=200, help="DPI of the page render the px boxes came from")
    p.add_argument("--from-recognize", help="JSON saved from `recognize` on the page render; uses its boxes")
    p.add_argument("--max-dist", type=float, help="max label distance from a structure, in points")
    p.set_defaults(func=cmd_annotate)

    p = sp.add_parser("export", help="write compound records to CSV/SDF/JSON")
    p.add_argument("records", help='JSON list: [{"id": "3a", "smiles": "...", "yield": "92", "ee": "95", ...}] or "-"')
    p.add_argument("-o", "--out", required=True, help="output file (.csv, .sdf or .json)")
    p.add_argument("--format", choices=["csv", "sdf", "json"])
    p.set_defaults(func=cmd_export)

    p = sp.add_parser("lookup", help="PubChem lookup (network)")
    p.add_argument("query", help="InChIKey, name, or SMILES (with --smiles)")
    p.add_argument("--smiles", action="store_true")
    p.set_defaults(func=cmd_lookup)

    p = sp.add_parser("engines", help="list available OCSR engines")
    p.set_defaults(func=cmd_engines)

    a = ap.parse_args()
    if a.cmd in RDKIT_COMMANDS and not HAVE_RDKIT:
        out({**RDKIT_MISSING, "command": a.cmd, "detail": RDKIT_ERROR})
        sys.exit(3)
    a.func(a)


if __name__ == "__main__":
    main()
