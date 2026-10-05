"""Tests for the commands added in 0.2.0: cdxml, docx-figures, grid, and the
RDKit-missing guard. Run: python3 -m pytest tests/"""
import json, os, subprocess, sys, tempfile, zipfile
HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "skills", "chem-vision", "scripts")
CV = os.path.join(SCRIPTS, "chemvision.py")
sys.path.insert(0, SCRIPTS)


def run(*args):
    r = subprocess.run([sys.executable, CV, *args], capture_output=True, text=True)
    return r.returncode, json.loads(r.stdout) if r.stdout.strip() else None


def _frag(fid, x, y):
    # ethanol C-C-O drawn left to right; an "O" atom label sits inside a node
    return f"""<fragment id="{fid}" BoundingBox="{x} {y} {x + 30} {y + 10}">
  <n id="{fid}1" p="{x} {y + 5}"/>
  <n id="{fid}2" p="{x + 15} {y}"/>
  <n id="{fid}3" p="{x + 30} {y + 5}" Element="8"><t p="{x + 28} {y + 8}" BoundingBox="{x + 27} {y + 1} {x + 33} {y + 9}"><s>O</s></t></n>
  <b id="{fid}4" B="{fid}1" E="{fid}2"/>
  <b id="{fid}5" B="{fid}2" E="{fid}3"/>
</fragment>"""


def _label(text, x, y):
    # ChemDraw style: bold compound number in its own run, line breaks kept inside the runs
    first, rest = text.split("\n", 1)
    runs = f'<s face="1">{first}\n</s><s>{rest}</s>'
    return f'<t p="{x} {y}" BoundingBox="{x} {y - 8} {x + 40} {y + 14}">{runs}</t>'


def test_cdxml_labels():
    doc = f"""<?xml version="1.0" encoding="UTF-8" ?>
<CDXML><page id="1">
{_frag(100, 50, 50)}{_label("3a" + chr(10) + "92% yield" + chr(10) + "95% ee", 50, 75)}
{_frag(200, 150, 50)}{_label("3b" + chr(10) + "88% (80%) yield" + chr(10) + "99% ee", 150, 75)}
</page></CDXML>"""
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "t.cdxml")
        open(p, "w").write(doc)
        code, res = run("cdxml", p)
        assert code == 0, res
        got = {s["compound_id"]: s for s in res["structures"]}
        assert set(got) == {"3a", "3b"}
        assert (got["3a"]["yield"], got["3a"]["ee"]) == ("92", "95")
        assert (got["3b"]["yield"], got["3b"]["ee"]) == ("88", "99")
        # atom labels inside nodes must not be read as captions
        assert "O" not in (res["unassigned_text"] or [])
        try:
            import rdkit  # noqa: F401
        except ImportError:
            return
        # With RDKit, each labelled structure should carry its SMILES, tied by fragment id.
        assert got["3a"].get("smiles") == "CCO", (
            "cdxml SMILES not mapped to fragments; RDKit may not set CDXML_FRAG_ID. notes: %s" % res["notes"])
        code, res = run("cdxml", p, "--id", "3b")
        assert [s["compound_id"] for s in res["structures"]] == ["3b"]


def test_docx_figures_png():
    from PIL import Image
    with tempfile.TemporaryDirectory() as d:
        png = os.path.join(d, "a.png")
        Image.new("RGB", (120, 90), "white").save(png)
        docx = os.path.join(d, "t.docx")
        W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
        R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
        A = "http://schemas.openxmlformats.org/drawingml/2006/main"
        with zipfile.ZipFile(docx, "w") as z:
            z.writestr("word/document.xml", f'<w:document xmlns:w="{W}" xmlns:r="{R}" xmlns:a="{A}"><w:body>'
                       f'<w:p><w:r><w:t>Scheme 2. Scope</w:t></w:r></w:p>'
                       f'<w:p><w:r><a:blip r:embed="rId9"/></w:r></w:p></w:body></w:document>')
            z.writestr("word/_rels/document.xml.rels",
                       '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                       '<Relationship Id="rId9" Target="media/image1.png"/></Relationships>')
            z.write(png, "word/media/image1.png")
        code, res = run("docx-figures", docx, "--outdir", os.path.join(d, "out"), "--grid")
        assert code == 0, res
        assert res["figure_count"] == 1
        f = res["figures"][0]
        assert os.path.exists(f["png"]) and os.path.exists(f["grid"]["grid_png"])
        assert "Scheme 2" in f["nearby_text"]


def test_grid():
    from PIL import Image
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "big.png")
        Image.new("RGB", (3000, 2000), "white").save(p)
        code, res = run("grid", p, "--max-width", "1000")
        assert code == 0
        assert res["full_size"] == [3000, 2000] and res["view_size"][0] == 1000
        assert os.path.exists(res["grid_png"])


def test_rdkit_missing_guard():
    """Simulate an environment without RDKit: validation commands must exit 3
    with a JSON explanation instead of crashing, and engines must still run."""
    code = ("import sys, runpy; sys.modules['rdkit'] = None; sys.argv = ['chemvision.py'] + sys.argv[1:]; "
            f"runpy.run_path({CV!r}, run_name='__main__')")
    r = subprocess.run([sys.executable, "-c", code, "check", "CCO"], capture_output=True, text=True)
    assert r.returncode == 3, r.stderr
    assert json.loads(r.stdout)["error"].startswith("RDKit is not available")
    r = subprocess.run([sys.executable, "-c", code, "engines"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    eng = json.loads(r.stdout)
    assert eng["rdkit"] is None and "rdkit_missing" in eng
