# chem-vision

A Claude Code plugin that lets Claude read chemical structures and reaction
schemes from images, figures and PDFs and turn them into SMILES, formulas and
plain-language descriptions.

## How it works

Claude reads the drawing with its own vision, guided by the `chem-vision`
skill. A helper CLI (`skills/chem-vision/scripts/chemvision.py`) covers the
parts vision alone is bad at:

| Command | What it does |
|---|---|
| `pages` | Renders PDF pages to PNG, extracts embedded figures and page text |
| `crop` | Crops and upscales a region so small structures can be read closely |
| `recognize` | Second opinion from OCSR: OSRA run as a 5-pass preprocessing ensemble with voting; MolScribe/DECIMER used if installed |
| `build` | Builds a molecule from atom positions + wedge/hash bonds read off the drawing, so RDKit assigns R/S |
| `check` | RDKit validation: canonical SMILES, formula, MW, InChIKey, stereocentres, functional groups, misread warnings |
| `render` | Draws candidates next to the original image for a visual check |
| `compare` | Agreement between candidates (InChIKey, skeleton, Tanimoto) |
| `reaction` | Parses a reaction SMILES, checks atom balance, records conditions/yield, draws the step |
| `annotate` | Compound IDs, yields, ee/er, dr, conversion, time, temperature next to structures; text-layer tables (Entry / Yield / ee columns) joined to structures by row |
| `export` | Writes compound records to CSV, SDF or JSON, flagging inconsistent values |
| `lookup` | PubChem name lookup by InChIKey/name/SMILES (needs network) |

## Install

In Claude Code:

```
/plugin marketplace add nathanohare/claude-chem-vision
/plugin install chem-vision@chem-vision
```

Then install the Python and OSRA dependencies once:

```bash
./setup.sh          # RDKit, PyMuPDF, Pillow, OSRA
./setup.sh --ml     # optional: MolScribe + DECIMER (large, download weights on first use)
```

Ask Claude things like "extract every compound with its yield and ee from
Table 2 of this PDF as a CSV" or "what is the product of step 3 in Scheme 1?".

## Test results (prototype, 2026-10-05)

`tests/make_samples.py` builds a synthetic set: 8 single structures in varied
styles (clean, black-and-white, low-res JPEG, abbreviations, rotated, noisy
scan, two complex natural products), a two-step reaction scheme with reagents
over the arrows, and a mock journal page PDF with a scheme and a figure.

A fresh Claude given only the skill and the images (no answer key) got
**18/18 structures exactly right including stereochemistry**, plus both
reaction steps with conditions and yields (`tests/blind_run/`, scored with
`tests/score_results.py`). OSRA alone got 3/8 single structures exactly
right (5/8 with the right skeleton) and nothing on the noisy scan
(`tests/run_osra_eval.py`).

Caveats: the images were drawn with RDKit, not ChemDraw scans from real
papers, and the molecules are well known, so expect lower accuracy on real
literature figures. Getting stereo right on paclitaxel and strychnine needed
the `build` step: from memory alone, the tester's first answers had wrong
configurations.

### Real papers (2026-10-05)

Four journal articles (Sonogashira coupling, Rh-catalysed asymmetric
1,4-additions, aminocyclopentadiene ring contractions, piperidine fragments;
not included in this repo) were each read end to end by a fresh Claude with
the skill: about 160 compound records with IDs, yields, ee and dr, plus
reaction SMILES for each scheme step. There is no answer key, but on the
vector-figure paper with a text layer, the automatic structure + label +
yield + ee join for Table 2 matched the hand-checked reading on all 5 rows.
Problems those runs exposed (labels like `3h` read as "3 h", symbol-font
glyphs, inverted CCITT figures, table columns, OSRA merge bugs) are fixed
and covered by `tests/test_annotations.py`.

## Known limits

- MolScribe and DECIMER adapters are written but untested here: the build
  environment could not reach Hugging Face or Zenodo for model weights.
- Markush structures, polymers and organometallics are described in words
  with partial SMILES.
- OSRA is unreliable on noisy scans and complex drawings; it is a second
  opinion, never the answer.
- Labels and values only come out of the PDF text layer for vector figures
  with real text. Many publishers outline figure text or embed raster
  images; there Claude reads the labels visually.
- Perspective drawings (chairs, bicyclics drawn in 3D) need Claude to reason
  about relative stereo; confidence is lower there.
