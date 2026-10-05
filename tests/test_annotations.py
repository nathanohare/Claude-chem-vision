"""Unit tests for label/value parsing and the PDF scope-table pipeline.
Run: python3 -m pytest tests/  (or python3 tests/test_annotations.py)"""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "skills", "chem-vision", "scripts")
sys.path.insert(0, SCRIPTS)
from annotations import annotate_text  # noqa: E402

CASES = [
    ("3a, 92%, 95% ee", {"ids": ["3a"], "yield": "92", "ee": "95"}),
    ("4b (85% yield, 97:3 er)", {"ids": ["4b"], "yield": "85", "er": "97:3", "ee_from_er": 94.0}),
    ("12c: 71%, >20:1 dr, 90% ee", {"ids": ["12c"], "yield": "71", "dr": ">20:1", "ee": "90"}),
    ("(±)-5, 64% (NMR yield)", {"ids": ["(±)-5"], "yield_nmr": "64"}),
    ("S7, 88%, E/Z = 10:1", {"ids": ["S7"], "yield": "88", "ez": "10:1"}),
    ("6aa 40% conv., 12 h, −20 °C", {"ids": ["6aa"], "conversion": "40", "time": "12 h", "temperature": "-20 °C"}),
    ("Conditions: 1 (0.5 mmol), Ru (1 mol%), HCO2H/Et3N, 30 °C, 24 h.",
     {"ids": [], "scale": "0.5 mmol", "loading": "1 mol%", "temperature": "30 °C", "time": "24 h"}),
    ("Table 1, entry 3: 3a, 92%", {"ids": ["3a"], "yield": "92"}),
    ("3d", {"ids": ["3d"]}),
]


def test_parse():
    for text, want in CASES:
        got = annotate_text(text)
        for k, v in want.items():
            assert got.get(k) == v, f"{text!r}: {k}={got.get(k)!r}, expected {v!r}"
        if "ids" in want:
            assert "time" not in got or "time" in want, f"{text!r}: spurious time {got.get('time')}"


def test_scope_pdf():
    pdf = os.path.join(HERE, "images", "scope_table.pdf")
    if not os.path.exists(pdf):
        subprocess.run([sys.executable, os.path.join(HERE, "make_scope_table.py")], cwd=HERE, check=True)
    truth = json.load(open(os.path.join(HERE, "scope_truth.json")))
    boxes = json.dumps([v["box_pt"] for v in truth.values()])
    cli = os.path.join(SCRIPTS, "chemvision.py")
    r = json.loads(subprocess.run([sys.executable, cli, "annotate", "--pdf", pdf, "--page", "1",
                                   "--boxes", boxes, "--units", "pt"],
                                  capture_output=True, text=True, check=True).stdout)
    got = {s["compound_id"]: s for s in r["structures"]}
    assert list(got) == list(truth)
    for cid, t in truth.items():
        joined = " ".join(t["labels"])
        for key in ("yield", "ee", "er", "dr"):
            v = got[cid].get(key)
            if v is not None:
                assert v in joined, f"{cid}: {key}={v} not in {joined!r}"
        assert got[cid].get("yield") is not None


if __name__ == "__main__":
    test_parse(); test_scope_pdf(); print("all passed")


def test_text_table():
    """Entry/Yield/ee table with footnote marks, '>' set in a symbol-like
    font, and a footnote line that must end the table."""
    import pymupdf as fitz
    from annotations import find_tables, group_lines, page_spans
    doc = fitz.open()
    page = doc.new_page(width=612, height=792)
    for x, h in [(50, "Entry"), (100, "Ligand"), (200, "Conv. (%)b"), (280, "Yield (%)c"), (360, "ee (%)")]:
        page.insert_text((x, 100), h, fontsize=8)
    rows = [("1", "L1", ">99", "92", "95"), ("2", "L2", "87", "63", "80"), ("3d", "L3", "40", "31", "12")]
    for i, (e, lig, conv, y, ee) in enumerate(rows):
        yy = 115 + 12 * i
        for x, v in [(50, e), (100, lig), (200, conv), (280, y), (360, ee)]:
            page.insert_text((x, yy), v, fontsize=8)
    page.insert_text((50, 160), "a Reaction conditions: substrate (0.2 mmol), 24 h.", fontsize=7)
    t = find_tables(group_lines(page_spans(page)))
    assert len(t) == 1
    got = [r["fields"] for r in t[0]["rows"]]
    assert [g["entry"] for g in got] == ["1", "2", "3d"]
    assert [g["yield"] for g in got] == ["92", "63", "31"]
    assert [g["ee"] for g in got] == ["95", "80", "12"]
    assert got[0]["conversion"] == ">99"
