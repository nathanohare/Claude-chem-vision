"""Score the OSRA engine alone against truth.json (exact InChIKey match + Tanimoto)."""
import json, subprocess, sys, os
here = os.path.dirname(os.path.abspath(__file__))
cli = os.path.join(here, "..", "skills", "chem-vision", "scripts", "chemvision.py")
truth = json.load(open(os.path.join(here, "truth.json")))
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
G = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
def fp(s): return G.GetFingerprint(Chem.MolFromSmiles(s))
for name in sorted(os.listdir(os.path.join(here, "images"))):
    key = os.path.splitext(name)[0]
    if key not in truth: continue
    r = json.loads(subprocess.run([sys.executable, cli, "recognize", os.path.join(here, "images", name), "--engine", "osra"], capture_output=True, text=True).stdout)
    t = Chem.MolToInchiKey(Chem.MolFromSmiles(truth[key]))
    best = None
    for h in r["results"]:
        if not h.get("valid"): continue
        sim = DataStructs.TanimotoSimilarity(fp(truth[key]), fp(h["canonical_smiles"]))
        k = h["inchikey"] or ""; exact = k == t; skel = k[:14] == t[:14]
        if best is None or (exact, skel, sim) > best[:3]: best = (exact, skel, sim, h["canonical_smiles"])
    if best: print(f"{key:18s} exact={best[0]!s:5s} skeleton={best[1]!s:5s} tanimoto={best[2]:.2f}  {best[3]}")
    else: print(f"{key:18s} no valid result  {r['notes']}")
