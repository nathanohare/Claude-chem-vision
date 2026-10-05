"""Score a results.json (from a blind run of the skill) against truth.json.
Usage: python3 score_results.py results.json blindmap.txt"""
import json, sys, re
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
G = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
truth = json.load(open(__file__.rsplit("/", 1)[0] + "/truth.json"))
res = json.load(open(sys.argv[1]))
amap = dict(re.match(r"(\S+) <- (\S+)", l).groups() for l in open(sys.argv[2]) if "<-" in l)
EXPECT = {"scheme_suzuki": ["scheme_1", "scheme_2", "scheme_3", "scheme_4"],
          "mock_paper": ["scheme_1", "scheme_2", "scheme_3", "scheme_4", "boc_ester", "testosterone"]}
def key(s):
    m = Chem.MolFromSmiles(s or "")
    return (Chem.MolToInchiKey(m), G.GetFingerprint(m)) if m else (None, None)
rows = []
for sample, orig in amap.items():
    stem = orig.rsplit(".", 1)[0]
    want = EXPECT.get(stem, [stem])
    got = [key(s.get("smiles")) for s in res.get(sample, {}).get("structures", [])]
    for w in want:
        tk, tf = key(truth[w])
        best = max(((k == tk, k is not None and k[:14] == tk[:14], DataStructs.TanimotoSimilarity(tf, f) if f else 0) for k, f in got), default=(False, False, 0))
        rows.append((orig, w, *best))
for r in rows: print(f"{r[0]:22s} {r[1]:14s} exact={r[2]!s:5s} skeleton={r[3]!s:5s} tanimoto={r[4]:.2f}")
n = len(rows); print(f"\nexact {sum(r[2] for r in rows)}/{n}, skeleton {sum(r[3] for r in rows)}/{n}, mean tanimoto {sum(r[4] for r in rows)/n:.2f}")
