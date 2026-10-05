"""Compound annotations: compound IDs, yields, ee/er, dr and similar values
that papers print next to structures (substrate-scope tables, schemes).

Two sources of text:
  * vector PDF text (exact, with positions and font flags) via PyMuPDF;
  * free text Claude transcribes from a raster figure.
Both go through the same parser, so results look the same either way.
"""
import re
import statistics

# Unicode variants papers use for minus, colon, etc.
_NORM = str.maketrans({"−": "-", "–": "-", "—": "-", "：": ":", "∶": ":",
                       " ": " ", " ": " ", " ": " ", "′": "'", "’": "'"})

NUM = r"\d{1,3}(?:\.\d+)?"
RATIO = rf"(?:[<>≥≤]\s*)?{NUM}\s*:\s*{NUM}"

# Order matters: specific patterns first, then the generic bare percentage.
# Units that collide with compound labels (3h, 8g, 1d) need a space before
# the unit: "24 h" is a time, "3h" is a compound.
PATTERNS = [
    ("ee", rf"(?P<v>[<>≥≤]?\s*{NUM})\s*%\s*ee\b"),
    ("ee", rf"\bee\s*(?:[=:]|of)?\s*(?P<v>[<>≥≤]?\s*{NUM})\s*%"),
    ("ee", rf"enantiomeric\s+excess(?:\s+of)?\s*(?P<v>[<>≥≤]?\s*{NUM})\s*%"),
    ("de", rf"(?P<v>[<>≥≤]?\s*{NUM})\s*%\s*de\b"),
    ("er", rf"(?P<v>{RATIO})\s*(?:e\.?r\.?|er)\b"),
    ("er", rf"\b(?:e\.?r\.?|er)\s*[=:]?\s*(?P<v>{RATIO})"),
    ("dr", rf"(?P<v>{RATIO})\s*(?:d\.?r\.?|dr)\b"),
    ("dr", rf"\b(?:d\.?r\.?|dr)\s*[=:]?\s*(?P<v>{RATIO})"),
    ("rr", rf"(?P<v>{RATIO})\s*(?:r\.?r\.?|rr)\b"),
    ("rr", rf"\b(?:r\.?r\.?|rr)\s*[=:]?\s*(?P<v>{RATIO})"),
    ("ez", rf"(?P<v>{RATIO})\s*(?:E\s*/\s*Z|Z\s*/\s*E)\b"),
    ("ez", rf"\b(?:E\s*/\s*Z|Z\s*/\s*E)\s*[=:]?\s*(?P<v>{RATIO})"),
    ("conversion", rf"(?P<v>{NUM})\s*%\s*conv(?:ersion|\.)?"),
    ("conversion", rf"\bconv(?:ersion|\.)?\s*(?:\(%\))?\s*[=:]?\s*(?P<v>{NUM})\s*%"),
    ("yield_nmr", rf"(?P<v>{NUM})\s*%\s*(?:\(\s*)?(?:NMR|1H NMR|GC|HPLC)\s*(?:yield)?\)?"),
    ("yield", rf"(?P<v>[<>≥≤]?\s*{NUM})\s*%\s*(?:isolated\s+)?yield"),
    ("yield", rf"\byield\s*[=:]?\s*(?P<v>{NUM})\s*%"),
    ("loading", r"(?P<v>\d+(?:\.\d+)?)\s*mol\s*%"),
    ("temperature", r"(?P<v>-?\d{1,3})\s*°\s*C\b"),
    ("temperature", r"\b(?P<v>rt|r\.t\.|room temperature|reflux)\b"),
    ("time", r"(?<![\w.-])(?P<v>\d+(?:\.\d+)?)\s*(?P<u>min|hr|hrs)\b"),
    ("time", r"(?<![\w.-])(?P<v>\d+(?:\.\d+)?)\s+(?P<u>h|days?)\b"),
    ("scale", r"(?<![\w.-])(?P<v>\d+(?:\.\d+)?)\s*(?P<u>mmol|mg)\b"),
    ("scale", r"(?<![\w.-])(?P<v>\d+(?:\.\d+)?)\s+(?P<u>mol|g)\b(?!\s*%)"),
]
_CASE_SENSITIVE = ("ez", "time", "scale", "temperature")
_COMPILED = [(k, re.compile(p, 0 if k in _CASE_SENSITIVE else re.I)) for k, p in PATTERNS]
_PCT = re.compile(rf"(?P<v>[<>≥≤]?\s*{NUM})\s*%")

# Footnote markers: "90%a", "8k[c]", "8k^c", "[a]" -> captured, then removed.
_FOOT = re.compile(r"(?<=%)(?P<f>[a-h])\b|\[(?P<g>[a-h])\]|\^(?P<h>[a-h])\b")

# Compound labels: 3, 3a, 12b', S1, 4aa, (±)-5, ent-7, (R)-3b, cis-5a, trans-5h
ID_RE = re.compile(r"^(?:\((?:[±+-]|[RS](?:,\s*[RS])*|rac)\)-|ent-|rac-|rel-|cis-|trans-|meso-|syn-|anti-|"
                   r"\([EZ]\)-)?(?:S|SI-)?\d{1,3}[a-z]{0,3}'*$")
_NOT_ID_BEFORE = re.compile(r"(?i)(table|scheme|figure|fig\.?|entry|entries|ref\.?|refs\.?|eq\.?|equiv\.?|"
                            r"step|chart|page|ca\.?|n\s*=)$")


def normalize(text):
    return text.translate(_NORM)


def split_footnotes(text):
    """Remove footnote markers, returning (clean_text, [markers])."""
    notes = []

    def grab(m):
        notes.append(m.group("f") or m.group("g") or m.group("h"))
        return ""

    return _FOOT.sub(grab, text), notes


def parse_values(text):
    """Extract labelled values from a snippet like '3a, 92%, 95:5 er'.
    A percentage with no qualifier is taken as the yield (the usual
    convention in scope tables) and flagged, but only in short label-like
    snippets: in running text a bare % is too ambiguous."""
    t = normalize(text)
    found, spans = {}, []

    def taken(a, b):
        return any(a < y and b > x for x, y in spans)

    # A parenthesised second result, e.g. "97% (58%: Pd = 3 mol %)", is kept
    # whole as an alternative instead of leaking its values into this entry.
    for m in re.finditer(r"\((?=[^()]*%)[^()]*\)", t):
        inner = m.group(0)[1:-1]
        if re.fullmatch(rf"\s*{NUM}\s*%\s*", inner) and not _PCT.search(t[:m.start()]):
            continue  # "(77%)" alone is the value itself
        if _PCT.search(t[:m.start()]):
            found.setdefault("alternative_results", []).append(inner.strip())
            spans.append(m.span())

    for key, rx in _COMPILED:
        for m in rx.finditer(t):
            if taken(*m.span()):
                continue
            v = re.sub(r"\s+", "", m.group("v"))
            if "u" in rx.groupindex and m.group("u"):
                v = f"{v} {m.group('u')}"
            elif key == "temperature" and v[-1].isdigit():
                v = f"{v} °C"
            elif key == "loading":
                v = f"{v} mol%"
            found.setdefault(key, v)
            spans.append(m.span())
    labelish = len(t.split()) <= 8
    for m in _PCT.finditer(t):
        if taken(*m.span()):
            continue
        v = re.sub(r"\s+", "", m.group("v"))
        if "yield" not in found and labelish:
            found["yield"] = v
            found["yield_assumed"] = True  # bare % taken as yield
        else:
            found.setdefault("other_percent", []).append(v)
        spans.append(m.span())
    # er -> ee convenience (for the major enantiomer).
    if "er" in found and "ee" not in found:
        try:
            a, b = [float(x) for x in re.sub(r"[<>≥≤]", "", found["er"]).split(":")]
            found["ee_from_er"] = round(abs(a - b) / (a + b) * 100, 1)
        except ValueError:
            pass
    return found


def find_ids(text):
    """Candidate compound labels in free text (tokens shaped like 3a, 12, S4)."""
    out = []
    toks = [t for t in re.split(r"[\s,;]+", normalize(text)) if t]
    for n, raw in enumerate(toks):
        if n and _NOT_ID_BEFORE.search(toks[n - 1]):
            continue
        if re.fullmatch(r"\(\d+\)", raw):  # "(4)": an equation number
            continue
        if re.fullmatch(r"\d+\.", raw):  # "2. H2O": a numbered step
            continue
        tok = raw.rstrip(".,;:")
        tok = tok.strip("()[]") if not (tok.startswith("(") and ")-" in tok) else tok
        if tok and ID_RE.match(tok) and not re.fullmatch(r"\d{4}", tok):
            out.append(tok)
    return out


def rank_ids(ids, bold=()):
    """Most label-like first: bold, then letter-suffixed (3c), then bare."""
    def score(i):
        return (i not in bold, not re.search(r"\d[a-z]", i), ids.index(i))
    return sorted(dict.fromkeys(ids), key=score)


def annotate_text(text):
    clean, notes = split_footnotes(text)
    vals = parse_values(clean)
    # IDs inside a parenthesised second result ("(40% recovered 8c)") are
    # not this entry's label.
    id_text = clean
    if "alternative_results" in vals:
        id_text = re.sub(r"\((?=[^()]*%)[^()]*\)", " ", normalize(clean))
    ids = find_ids(id_text)
    # Drop tokens that are really numbers inside values ("95" of "95:5 er").
    value_nums = set(re.findall(r"\d+(?:\.\d+)?", " ".join(
        str(v) for k, v in vals.items() if k not in ("yield_assumed",))))
    ids = [i for i in ids if i not in value_nums]
    if len(clean.split()) > 8:  # running text: bare numbers are rarely labels
        ids = [i for i in ids if not i.isdigit()]
    rec = {"text": text, "ids": rank_ids(ids), **vals}
    if notes:
        rec["footnotes"] = notes
    return rec


# ------------------------------------------------------------- PDF geometry

# Elsevier and some other publishers set symbols in "Math Pi" fonts whose
# glyph codes look like ordinary letters.
_MATHPI = {"8": "°", "O": ">", "!": "<", "Z": "=", "K": "-", "C": "+", "G": "±", "P": "×",
           "R": "≥", "%": "≤"}
_BOLD_FONT = re.compile(r"(?i)(bold|black|heavy|semibold|demi|-b$|-bd|-bi$|-bold|,bold)")


def _is_mathpi(font):
    return bool(re.search(r"(?i)(MPi|MathPi|MathematicalPi|UniversalMath)", font))


def page_spans(page):
    """Text spans with bbox (PDF points), bold flag, and a flag for sub- and
    superscripts (smaller than the surrounding line)."""
    spans = []
    d = page.get_text("dict")
    all_sizes = [s["size"] for b in d["blocks"] for l in b.get("lines", []) for s in l["spans"]
                 if s["text"].strip()]
    body = statistics.median(all_sizes) if all_sizes else 0
    for b in d["blocks"]:
        for line in b.get("lines", []):
            sizes = [s["size"] for s in line["spans"] if s["text"].strip()]
            big = max(sizes + [body]) if sizes else 0
            prev_font = None
            raw = line["spans"]
            for k, s in enumerate(raw):
                txt = s["text"].strip()
                if not txt:
                    continue
                if _is_mathpi(s["font"]) and len(txt) <= 2:
                    txt = "".join(_MATHPI.get(c, c) for c in txt)
                elif txt == "8" and s["font"] != prev_font and k + 1 < len(raw) \
                        and raw[k + 1]["text"].strip().startswith("C"):
                    txt = "°"  # degree sign set in a symbol font
                elif txt in ("O", "!", "Z", "G") and k + 1 < len(raw) \
                        and raw[k + 1]["font"] != s["font"] \
                        and re.match(r"\s*\d", raw[k + 1]["text"]) \
                        and abs(raw[k + 1]["size"] - s["size"]) < 0.5:
                    txt = _MATHPI[txt]  # ">99" with ">" from a symbol font
                prev_font = s["font"]
                bold = bool(s["flags"] & 16) or bool(_BOLD_FONT.search(s["font"]))
                script = bool(s["flags"] & 1) or (big and s["size"] < 0.8 * big)
                spans.append({"text": txt, "bbox": [round(v, 1) for v in s["bbox"]],
                              "bold": bold, "script": script, "size": round(s["size"], 1)})
    return spans


def group_lines(spans, ygap=2.5):
    """Merge spans sitting on the same baseline into lines (ChemDraw text
    often arrives as one span per run: '3a', ', ', '92%', ...). Sub/superscripts
    join the preceding text without a space (BF3K, 90%a) and are never taken
    as labels on their own. A large horizontal gap starts a new line, so the
    two columns of a page stay apart."""
    body = statistics.median([s["size"] for s in spans]) if spans else 10
    spans = sorted(spans, key=lambda s: (round(s["bbox"][1] / 3), s["bbox"][0]))
    lines = []
    for s in spans:
        for L in reversed(lines):
            same_row = abs(L["bbox"][1] - s["bbox"][1]) < ygap + (4 if s["script"] else 0)
            gap = s["bbox"][0] - L["bbox"][2]
            if same_row and s["script"] and L["bbox"][0] < s["bbox"][0] < L["bbox"][2]:
                # subscript inside an already-merged run ("K PO4" + "3"):
                # splice it in at the matching character position
                frac = (s["bbox"][0] - L["bbox"][0]) / max(1e-6, L["bbox"][2] - L["bbox"][0])
                pos = round(frac * len(L["text"]))
                while 0 < pos < len(L["text"]) and L["text"][pos - 1] == " ":
                    pos -= 1
                rest = L["text"][pos:].lstrip(" ")
                L["text"] = L["text"][:pos] + s["text"] + rest
                break
            if same_row and -1 <= gap < 1.5 * max(s["size"], body * 0.8):
                glue = "" if s["script"] or L["text"].endswith(("(", "-")) else " "
                L["text"] += glue + s["text"]
                L["bbox"][2] = max(L["bbox"][2], s["bbox"][2])
                L["bbox"][3] = max(L["bbox"][3], s["bbox"][3])
                if s["bold"] and not s["script"]:
                    L["bold_parts"].append(s["text"])
                break
        else:
            lines.append({"text": s["text"], "bbox": list(s["bbox"]), "size": s["size"],
                          "bold_parts": [s["text"]] if s["bold"] and not s["script"] else []})
    return lines


# ------------------------------------------------------------------ tables

_HEADER_KW = re.compile(r"(?i)^(entry|run|yield|ee|er|dr|d\.r\.|e\.r\.|conv|time|t\s*\(|temp|solvent|"
                        r"catalyst|cat\.|ligand|base|additive|substrate|product|selectivity|ratio|"
                        r"aryl|alkyne|alkene|amine|ketone|aldehyde|.*\(%\)|.*mol\s*%|[A-Z]?\d?\s*/\s*\d)")
_FIELD = [("entry", r"(?i)^(entry|run)"), ("product_id", r"(?i)^product"),
          ("substrate_id", r"(?i)^(aryl|substrate|starting|alkyne|alkene|amine|ketone|aldehyde|reactant|sm\b)"), ("yield", r"(?i)yield"), ("ee", r"(?i)^ee[a-h]?\b|\bee[a-h]?\s*\(|e\.e\."),
          ("er", r"(?i)^er[a-h]?\b|\ber[a-h]?\s*\(|e\.r\."), ("dr", r"(?i)^dr[a-h]?\b|\bdr[a-h]?\s*\(|d\.r\."),
          ("conversion", r"(?i)conv"), ("time", r"(?i)^time|^t\s*\((h|min)"),
          ("temperature", r"(?i)^temp|^T\s*\(°"), ("loading", r"(?i)mol\s*%")]


def header_field(h):
    h = re.sub(r"(?<=\))[a-h]$|(?<=%)[a-h]$", "", h.strip())  # footnote marks
    for field, rx in _FIELD:
        if re.search(rx, h):
            return field
    return h


def _rows_by_y(lines, tol=3.0):
    rows = []
    for L in sorted(lines, key=lambda L: L["bbox"][1]):
        cy = (L["bbox"][1] + L["bbox"][3]) / 2
        if rows and abs(cy - rows[-1]["cy"]) < tol:
            rows[-1]["lines"].append(L)
        else:
            rows.append({"cy": cy, "lines": [L]})
    return rows


def find_tables(lines):
    """Text-layer tables: a header row with at least three column titles, at
    least one of them a results/condition keyword, then body rows anchored on
    the first column (usually Entry). Returns rows with raw cells and
    parsed fields (yield, ee, conversion, ...)."""
    tables = []
    ys = _rows_by_y(lines)
    for i, hdr in enumerate(ys):
        cells = sorted((L for L in hdr["lines"] if not re.fullmatch(r"[a-h]", L["text"].strip())),
                       key=lambda L: L["bbox"][0])  # lone footnote letters are not columns
        if len(cells) < 3 or sum(bool(_HEADER_KW.match(c["text"])) for c in cells) < 2:
            continue
        if not any(header_field(c["text"]) in ("yield", "ee", "er", "dr", "conversion", "entry") for c in cells):
            continue
        x_starts = [c["bbox"][0] for c in cells]
        x_end = max(c["bbox"][2] for c in cells) + 40
        names = [c["text"] for c in cells]

        def col_of(L):
            cx = L["bbox"][0]
            j = max([k for k, x in enumerate(x_starts) if x <= cx + 4] or [0])
            return j

        body, last_y = [], hdr["cy"]
        for row in ys[i + 1:]:
            if row["cy"] - last_y > 45:
                break
            wide = [L for L in row["lines"] if L["bbox"][2] - L["bbox"][0] > 0.6 * (x_end - x_starts[0])]
            if wide or any(re.match(r"(?i)^(table|scheme|figure)\b", L["text"]) for L in row["lines"]):
                break
            if any(_rows_by_y([L])[0] and sum(bool(_HEADER_KW.match(c["text"])) for c in row["lines"]) >= 3
                   for L in row["lines"][:1]):
                break  # next table's header
            body.append(row)
            last_y = row["cy"]
        # Rows start where the first column has a value. When that column is
        # "Entry", a value that does not look like an entry number (a footnote,
        # body text) ends the table.
        entry_col = header_field(names[0]) == "entry"
        kept = []
        for row in body:
            first = [L for L in row["lines"] if col_of(L) == 0 and x_starts[0] - 10 <= L["bbox"][0]]
            if entry_col and first and not re.fullmatch(r"\d{1,3}[a-h]?", first[0]["text"].strip()):
                break
            kept.append(row)
        cells_in = [L for row in kept for L in row["lines"]
                    if x_starts[0] - 10 <= L["bbox"][0] <= x_end]
        anchors = sorted((L for L in cells_in if col_of(L) == 0), key=lambda L: L["bbox"][1])
        out_rows = [{"cells": {names[0]: L["text"]}, "bbox": list(L["bbox"]), "lines": [(L, names[0])],
                     "cy": (L["bbox"][1] + L["bbox"][3]) / 2} for L in anchors]
        # Other cells go to the nearest row anchor (labels drawn under a tall
        # structure still belong to the row whose entry number is level with it).
        for L in sorted(cells_in, key=lambda L: (L["bbox"][1], L["bbox"][0])):
            j = col_of(L)
            if j == 0 or not out_rows:
                continue
            cy = (L["bbox"][1] + L["bbox"][3]) / 2
            r = min(out_rows, key=lambda r: abs(r["cy"] - cy))
            key = names[j]
            r["cells"][key] = (r["cells"][key] + " " + L["text"]) if key in r["cells"] else L["text"]
            r["lines"].append((L, key))
            r["bbox"] = [min(r["bbox"][0], L["bbox"][0]), min(r["bbox"][1], L["bbox"][1]),
                         max(r["bbox"][2], L["bbox"][2]), max(r["bbox"][3], L["bbox"][3])]
        for r in out_rows:
            r.pop("cy")
        for r in out_rows:
            parsed = {}
            for k, v in r["cells"].items():
                f = header_field(k)
                clean, notes = split_footnotes(v)
                vv = clean.strip()
                if f in ("yield", "ee", "conversion"):
                    m = re.search(r"[<>≥≤]?\s*\d+(?:\.\d+)?", vv)
                    vv = re.sub(r"\s+", "", m.group(0)) if m else vv
                if f in ("product_id", "substrate_id"):
                    ids = rank_ids(find_ids(vv))
                    vv = ids[0] if ids else vv
                parsed[f] = vv
                if notes:
                    parsed.setdefault("footnotes", []).extend(notes)
            r["fields"] = parsed
        if out_rows:
            tables.append({"columns": names, "header_lines": cells, "header_bbox": [x_starts[0], hdr["lines"][0]["bbox"][1],
                                                            x_end, hdr["lines"][0]["bbox"][3]],
                           "rows": out_rows})
    return tables


def _dist(box, line):
    """Distance from a structure box to a text line, favouring text below
    or beside the structure (where scope tables put labels)."""
    x0, y0, x1, y1 = box
    lx0, ly0, lx1, ly1 = line["bbox"]
    cx, cy = (lx0 + lx1) / 2, (ly0 + ly1) / 2
    dx = 0 if x0 <= cx <= x1 else min(abs(cx - x0), abs(cx - x1))
    dy = 0 if y0 <= cy <= y1 else min(abs(cy - y0), abs(cy - y1))
    if cy < y0:
        dy *= 1.6  # text above a structure is less often its label
    return (dx * dx + dy * dy) ** 0.5


def _line_records(lines):
    recs = []
    for L in lines:
        a = annotate_text(L["text"])
        bold_ids = [i for p in L["bold_parts"] for i in find_ids(split_footnotes(p)[0])]
        if not (a["ids"] or set(a) - {"text", "ids"}):
            continue
        a.update({"bbox": L["bbox"], "bold_ids": bold_ids})
        recs.append(a)
    return recs


def annotate_page(page, boxes=None, max_dist=None):
    """Find labels/values on a PDF page. If structure boxes (PDF points) are
    given, attach each text line to its nearest box within max_dist, then
    use table rows for what is left."""
    lines = group_lines(page_spans(page))
    tables = find_tables(lines)
    # Table cells are consumed by the table, except compound-label columns:
    # with structure boxes known, a label is tied to its drawing by position,
    # which is more reliable than guessing its row.
    in_table = set()
    for t in tables:
        in_table.update(id(L) for L in t.pop("header_lines"))
        for r in t["rows"]:
            for L, col in r.pop("lines"):
                if not (boxes and header_field(col) in ("product_id", "substrate_id")):
                    in_table.add(id(L))
    recs = _line_records([L for L in lines if id(L) not in in_table])
    if not boxes:
        return {"lines": recs, "tables": tables}
    # Was there any text near the structures at all? If not, the figure is
    # probably a raster image and the page text is just body text: attach
    # nothing rather than something wrong.
    near = [r for r in recs if min(_dist(b, r) for b in boxes) < 20]
    if not near and not tables:
        return {"structures": [{"box": b, "compound_id": None} for b in boxes], "unassigned_lines": [],
                "warning": ("No text-layer labels near the structures: the figure is probably a raster "
                            "image (or its labels are outlined paths). Read labels visually and use "
                            "--labels or --pair.")}
    res = assign(boxes, recs, max_dist)
    # Table rows level with structures: the row's values go to the rightmost
    # structure in the row band (the product, in reactant -> product tables).
    for t in tables:
        for r in t["rows"]:
            y0, y1 = r["bbox"][1], r["bbox"][3]
            cy = (y0 + y1) / 2
            band = [s for s in res["structures"] if s["box"][1] - 4 <= cy <= s["box"][3] + 4]
            if not band:
                continue
            band.sort(key=lambda s: s["box"][0])
            f = r["fields"]
            for s in band[:-1]:
                s.setdefault("table_entry", f.get("entry"))
                if f.get("substrate_id") and ID_RE.match(f["substrate_id"]) and not s.get("compound_id"):
                    s["compound_id"] = f["substrate_id"]
            prod = band[-1]
            for k, v in f.items():
                if k == "product_id":
                    if ID_RE.match(v) and not prod.get("compound_id"):
                        prod["compound_id"] = v
                    continue
                if k == "substrate_id":
                    continue
                prod.setdefault("table_entry" if k == "entry" else k, v)
            prod["assigned_by_table_row"] = True
            r["matched_structure"] = True
    res["tables"] = tables
    return res


def annotate_labels(boxes, labels, max_dist=None):
    """Same as annotate_page, for labels Claude read visually from a raster
    figure: [{"text": "3a, 92%, 95% ee", "box": [x0, y0, x1, y1]}, ...]
    with boxes in the same units as the structure boxes."""
    recs = []
    for lab in labels:
        a = annotate_text(lab["text"])
        a.update({"bbox": list(lab["box"]), "bold_ids": []})
        recs.append(a)
    return assign(boxes, recs, max_dist)


def assign(boxes, recs, max_dist=None):
    structs = [{"box": b, "annotations": []} for b in boxes]
    leftovers = []
    for r in recs:
        dists = [_dist(b, r) for b in boxes]
        j = min(range(len(boxes)), key=dists.__getitem__)
        h = boxes[j][3] - boxes[j][1]
        limit = max_dist if max_dist is not None else 0.6 * h + 12
        if dists[j] <= limit:
            r["assigned_by"] = "proximity"
            structs[j]["annotations"].append(r)
        else:
            leftovers.append(r)
    # Table rows: a value printed in its own column, level with a structure,
    # belongs to the nearest structure to its left in that row band.
    unassigned = []
    for r in leftovers:
        lx0, ly0, lx1, ly1 = r["bbox"]
        cy = (ly0 + ly1) / 2
        row = [j for j, b in enumerate(boxes) if b[1] <= cy <= b[3] and b[2] <= lx0 + 2]
        if row and set(r) - {"text", "ids", "bbox", "bold_ids", "footnotes"}:
            j = max(row, key=lambda j: boxes[j][2])
            r["assigned_by"] = "table_row"
            structs[j]["annotations"].append(r)
        else:
            unassigned.append(r)
    for s in structs:
        merged = {}
        for r in sorted(s["annotations"], key=lambda r: (r["assigned_by"] != "proximity",
                                                          r["bbox"][1], r["bbox"][0])):
            for k, v in r.items():
                if k in ("text", "bbox", "ids", "bold_ids", "assigned_by"):
                    continue
                if k == "footnotes":
                    merged.setdefault(k, []).extend(v)
                else:
                    merged.setdefault(k, v)
        bold = [i for r in s["annotations"] for i in r["bold_ids"]]
        ids = rank_ids([i for r in s["annotations"] for i in (r["bold_ids"] + r["ids"])], bold)
        s["compound_id"] = ids[0] if ids else None
        if len(ids) > 1:
            s["other_ids"] = ids[1:]
        s.update(merged)
        if any(r["assigned_by"] == "table_row" for r in s["annotations"]):
            s["assigned_by_table_row"] = True
        s["source_text"] = " | ".join(r["text"] for r in s["annotations"])
        del s["annotations"]
    return {"structures": structs, "unassigned_lines": unassigned}
