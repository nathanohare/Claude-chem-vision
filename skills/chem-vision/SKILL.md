---
name: chem-vision
description: Read chemical structures and reaction schemes from images, figures, PDFs, Word documents, ChemDraw files and scientific papers and turn them into SMILES with formulas, names, compound IDs, yields, ee/er and dr. Use whenever a user shares or points to a molecule drawing, skeletal formula, reaction scheme, synthetic route, or a paper/PDF page containing them and wants it identified, transcribed, checked, compared or explained.
---

# Reading chemical structures and schemes

You interpret the drawing with your own vision. The helper CLI supplies what
vision alone is bad at: zooming, an independent OCSR reading (OSRA, plus
MolScribe/DECIMER if installed), chemistry validation (RDKit), and re-drawing
your answer so you can check it against the original picture.

Perfect accuracy is not expected, especially for large natural products. Be
useful and honest: give your best structure, say how confident you are, and
point at the specific parts you are unsure about.

```
CV="python3 ${CLAUDE_PLUGIN_ROOT}/skills/chem-vision/scripts/chemvision.py"
```
(If `CLAUDE_PLUGIN_ROOT` is unset, use the path of this skill's `scripts/`
directory.) Every subcommand prints JSON. Run `$CV engines` once to see the RDKit
version and which recognizers are installed. If RDKit or PyMuPDF is missing,
run the plugin's `setup.sh` once; then follow **If RDKit is missing** below
for anything it could not install.

## If RDKit is missing

`engines` reports `"rdkit": null` with an `rdkit_missing` block when RDKit
could not be installed. If `setup.sh` or pip says "No matching distribution"
or returns 403, the package index is blocked by the environment's network
allowlist: do not retry with other install commands, mirrors or downloads.
Tell the user once, in one sentence, then carry on:

- Still available: `pages` (falls back to Poppler without PyMuPDF), `grid`,
  `crop`, `docx-figures`, `cdxml` (exact labels and values, no SMILES),
  `annotate --text/--labels/--pair`, `lookup`.
- Unavailable (they exit with code 3 and a JSON explanation): `recognize`,
  `check`, `render`, `build`, `compare`, `reaction`, `export`.
- Read structures visually and write SMILES by hand. Say in the report that
  the SMILES, formula and any R/S assignment are unvalidated, and give the
  stereodescriptor only with a stated confidence.
- The fix is to run the plugin where `pip install rdkit` succeeds (for
  example the user's own computer), or for an org admin to allow PyPI.

## Step 0: choose the best source

Work from the most exact source available, in this order:

1. **ChemDraw source (`.cdxml`)**: if the user has one, or one sits next to
   the document (same folder, a name like "schemes" or "scope"), use it:
   `$CV cdxml schemes.cdxml --id 3m` (or no `--id` for every labelled
   structure). Labels, yields and ee come back exact and tied to each
   structure by position; with RDKit, each structure's SMILES too. One
   canvas can hold several schemes or old versions, so the same ID may
   appear twice: the `notes` say so, and `box_pt` tells you which is which.
   Confirm with the user which scheme is current if it matters. `.cdx`
   (binary) is not supported by the XML reader; ask for a `.cdxml` export.
2. **Word document (`.docx`)**: `$CV docx-figures paper.docx --grid`
   extracts every embedded figure in document order, converts EMF/WMF
   (pasted ChemDraw) to PNG, trims page margins, and returns each figure's
   `title` ("Scheme 2. Scope of ...") and its real `text` when the figure
   has a text layer (EMF from ChemDraw usually does). Search `text` for the
   compound ID to get exact values. Do not convert the whole .docx to PDF:
   that loses figure text and resolution.
3. **PDF or image**: the page-rendering workflow below.

## Workflow

1. **Get images.**
   - PDF: `$CV pages paper.pdf --pages 3-5 --dpi 200 --text --outdir pages/`
     renders pages and returns page text (useful for compound numbers,
     captions and names the authors give). Add `--images` to also pull
     embedded raster figures at native resolution; when a figure is embedded,
     prefer that image over the page render (page renders can clip or
     overlap figures and are lower resolution).
   - Look at each page/figure image yourself (Read the PNG) and locate the
     structures, schemes and compound labels. Read the captions and
     footnotes in the page text too: they say what the numbers mean, give
     names for compounds that are not drawn, and list general conditions.
   - Find out early whether the figure labels are real text:
     `$CV annotate --pdf paper.pdf --page 4` lists every label-like line and
     any text-layer tables (Entry / Yield / ee columns) it finds. Journal
     figures come in three kinds: vector with real text (labels and values
     come out exactly), vector with outlined text (looks sharp but `annotate`
     finds nothing in the figure), and raster images (`pages --images`
     returns them). Only the first kind gives text-layer labels.

2. **Zoom in.** Small or dense structures are where reading errors come from.
   Find the structure on a grid first: `pages --grid` (or `docx-figures
   --grid`, or `$CV grid image.png`) writes a downscaled view with gridlines
   labelled `fraction|full-res-px`. **Choose crop boxes from that view and
   pass them as 0-1 fractions**: fractions are the same on every scale, so
   the first crop lands. Never estimate pixel coordinates from a downscaled
   view and apply them to the full-resolution image.
   Crop each structure (or each step of a scheme) and upscale:
   `$CV crop page.png --box 0.1,0.35,0.55,0.6 --scale 2 -o s1.png`
   (`--box` takes pixels or 0-1 fractions; `--pad` adds 10 px on each side
   by default, `--pad 0` for an exact box). Read the crop.

3. **Transcribe it yourself.** Write a SMILES by walking the drawing
   systematically: count ring atoms, note every heteroatom label, every double
   bond, every wedge/hash (stereo), charges, and abbreviations (Ph, Bn, Boc,
   Ts, OMe, TBS, Ac, Bz, Cbz, Fmoc, PMB, etc.; expand them to atoms). Implicit
   carbons sit at every vertex and line end. For R-groups/Markush structures,
   use `*` for attachment points and describe the R-group table in words.

   **Stereo: build from the drawing, not from memory.** For anything with
   wedges/hashes beyond one or two centres (steroids, sugars, natural
   products), do not hand-write @/@@. Place atoms at their approximate image
   coordinates and mark wedge/hash bonds, and let RDKit derive R/S:
   ```
   $CV build spec.json -o s1.mol
   # spec.json: {"atoms": [["C", 220, 22], ["N", 200, 185], ["O", 95, 185, -1], ...],
   #             "bonds": [[1, 2, 1], [2, 3, 2], [8, 9, 1, "wedge"], [12, 26, 1, "hash"], ...]}
   ```
   Atom indices are 1-based; a wedge/hash starts (narrow end) at the first
   atom of the bond; a 4th atom field is a formal charge; include explicit
   `["H", x, y]` atoms where the drawing shows wedged H. Coordinates only need
   to be roughly right (relative positions matter, scale does not). The
   saved `.mol` keeps the drawing's layout, so `render s1.mol --original ...`
   looks like the source and wedges can be compared one to one. Even if you
   recognise a famous molecule, check its stereo this way: memory is often
   wrong on configuration.

4. **Get a second opinion.** `$CV recognize s1.png` runs the OCSR engines
   (20-60 s per image or page; give the command a few minutes, or crop
   dense tables row by row, which is also more accurate).
   OSRA is run as an ensemble of preprocessing passes; `votes` says how many
   passes agreed. When an image holds several structures, results cover all
   of them: use `box_px` to tell which drawing each belongs to, or crop each
   structure first. Results containing `*` (an unread label, a compound number
   sitting on the drawing, or a cut-off edge) are listed last with a warning.
   OSRA is good on clean, simple drawings and often wrong on complex,
   rotated, noisy or abbreviation-heavy ones (typical failures: aryl groups
   read as cyclohexyl, alkynes lost, "Hex" read as a ring; results that look
   like this carry a `warning`); never trust it over a
   careful reading, but take a disagreement as a cue to look again. Options:
   `--adaptive`/`--jaggy` for scans, `--rotate DEG` for tilted drawings.

5. **Validate.** `$CV check "<smiles>"` (several at once is fine) gives
   validity, canonical SMILES, formula, MW, InChIKey, `stereocenters` (1-based
   atom numbers, like `build`; and `unassigned_stereocenters`), functional groups and `warnings` (radicals, odd
   elements, stray fragments). Fix invalid SMILES before going on. If the
   paper states a formula, mass or HRMS value, compare it: a formula mismatch
   means a misread.

6. **Check visually.** `$CV render "<your smiles>" "<osra smiles>" --original s1.png -o cmp.png`
   draws the original next to each candidate (SMILES get a fresh RDKit
   layout, which may be rotated or mirrored relative to the source; pass the
   `.mol` from `build` to keep the source layout). Read `cmp.png` and compare ring
   by ring, substituent by substituent, wedge by wedge. Correct and repeat
   until the drawing matches (at most a few rounds). `$CV compare a b`
   quantifies agreement (identical InChIKey, same skeleton ignoring stereo,
   `enantiomers`, Tanimoto). `--cols` sets panels per row for many candidates.

   **Racemic and relative stereo.** Wedges in a racemic or "rel" compound
   show relative configuration only. Draw it as shown (one enantiomer), say
   "racemic" or "relative configuration" in `notes`, and treat an
   `enantiomers: true` comparison as a match. Paired diastereomers (cis-5a /
   trans-5a) are separate records with the prefix kept in the ID. If a paper
   uses one label for two drawings (e.g. both enantiomers), keep the label and
   add a suffix (`3aa`, `3aa_ent`) with a note.

7. **Reaction schemes.** Identify reactants, products and the reagents and
   conditions over and under each arrow (read these as text: catalysts,
   solvents, temperature, time, yield). Build a reaction SMILES
   `reactants>agents>products` per step and run
   `$CV reaction "<rxn>" --label "1 + 2 -> 3" --conditions "Pd(PPh3)4, K2CO3, dioxane/H2O, 90 °C" --yield 87% -o step1.png`.
   `product_minus_reactant_atoms` shows what the arrow implies; check that it makes chemical sense for the stated
   reagents (e.g. a reduction should add H2; a Suzuki coupling loses B(OH)2
   and Br). Multi-step schemes: one reaction SMILES per arrow, keep compound
   numbers from the paper. Put the main product first; co-products or
   isomers after it are ignored in the atom difference.

8. **Compound IDs, yields, ee and other values.** Papers print a label and
   results next to each structure (`3a, 92%, 95% ee`, `4b (85% yield, 97:3
   er)`, `>20:1 dr`). Always capture them and tie each one to its structure.
   - PDF with a text layer (most journal PDFs, including vector ChemDraw
     figures): render the page, run `recognize` on it and save the JSON, then
     `$CV annotate --pdf paper.pdf --page 4 --from-recognize rec.json --dpi 200`
     (`--dpi` must match the page render `recognize` ran on). If you ran
     `recognize` on an embedded image from `pages --images` instead, pass
     `--image-rect` with that image's `page_rect_pt` and drop `--dpi`.
     Labels are tied to structures by position, and values printed in
     table columns (Entry, Yield (%), ee (%), Conv. (%), Pd (mol%) ...) are
     tied by row: the row's values go to the rightmost structure in that row
     (the product), and earlier structures in the row get its `table_entry`.
     Each structure comes back with `compound_id` (bold labels preferred),
     `yield`, `ee`, `er` (plus `ee_from_er`), `dr`, `de`, `rr`, `ez`,
     `conversion`, `yield_nmr`, `time`, `temperature`, `scale`, `loading` and
     the `source_text` it used. Lines it could not place are listed in
     `unassigned_lines` (often the table's general conditions), and the
     parsed `tables` are included. If the result has a `warning` that no
     labels are near the structures, the figure's text is not real text:
     switch to the raster route below. If OSRA misses
     or misreads a structure, pass your own boxes instead:
     `--boxes '[[x0,y0,x1,y1], ...]' --units px|pt|frac`.
   - Raster figures and scans (`annotate` says there is no text layer): read
     the labels yourself. Either give them with their approximate positions,
     in the same pixel space as the structure boxes, and let `annotate`
     attach them to the nearest structure:
     `--labels '[{"text": "3a, 92%, 95% ee", "box": [300,640,470,700]}, ...]'`
     or, simplest, pair each final structure with its label text directly:
     `$CV annotate --pair 3a.mol "3a, 92%, 95% ee" --pair "CCO..." "3b, 88%, 97:3 er" > compounds.json`
     (a structure is a SMILES or a `.mol` from `build`). The output is ready
     for `export`. `annotate --text "..."` parses text alone.
   - Tables with no structures (optimisation tables): `annotate --pdf`
     without boxes returns them as rows; make one record per entry with the
     product's SMILES and `entry` in the ID or notes.
   - Compounds named but not drawn (starting materials, "+ trans-5e, 22%"):
     include them from the text or another figure, with `confidence: medium`
     and a note saying where the structure came from.
   - A bare percentage next to a structure is taken as the isolated yield
     (`yield_assumed: true`); check the table footnotes for what the numbers
     mean (NMR yield, conversion, ee by HPLC/SFC, values in parentheses for a
     second entry) and say so. Footnote markers (`90%a`, `8k[c]`) are
     returned in `footnotes`; look them up and put what they say in `notes`.
     A parenthesised second result ("97% (58%: Pd = 3 mol%)") is kept in
     `alternative_results` rather than mixed into the entry.
   - Keep values as printed: an er stays in `er` (the tool adds
     `ee_from_er`); do not move it into `ee`.
   - General conditions (table header or footnote): put the shared text in a
     `conditions` field on every record and fill `temperature`, `time`,
     `scale` from it; when an entry prints its own value ("48 h"), that
     entry's value wins.
   - Your own reading of the structure always wins over OSRA's; the label and
     values from the text layer are exact, so prefer them over reading text
     visually.

9. **Export.** Collect one record per compound and write a table:
   `$CV export compounds.json -o scope.csv` (or `.sdf`, `.json`). Records are
   `{"id": "3a", "smiles": "...", "molfile": "3a.mol", "yield": "92",
   "ee": "95", "page": 4, "figure": "Table 2", "notes": "..."}`; any extra
   keys become columns. Export recomputes formula, MW and InChIKey from the
   SMILES, uses `molfile` (when given) for the SDF so it keeps the paper's
   layout, and lists `problems`: invalid SMILES, a dr on a structure with
   fewer than two stereocentres, ee on an achiral structure. Each problem
   means the structure or the value assignment needs another look. The JSON
   from `annotate --from-recognize` can be exported directly, but fix any
   structure you corrected first.

10. **Optional identification.** `$CV lookup <InChIKey>` asks PubChem for a
   name and IUPAC name (needs network). Only claim a common name if the
   lookup or the paper confirms it, or you are sure.

## What to report

For each structure, keep it short:
- the label used in the document (e.g. "3a"), SMILES, formula, MW, and the
  values printed with it (yield, ee/er, dr, ...);
- what it is (compound class, key functional groups, name if known);
- confidence: high / medium / low, and *where* the doubt is ("the C7
  stereocentre is hashed in the figure but blurry"; "OSRA disagrees on ring
  fusion"). Mention when OCSR and your reading agreed exactly; that is good
  evidence;
- for schemes: a step table (reactants, reagents/conditions, product, yield)
  plus the reaction SMILES, and a one-line description of each transformation.

When there are more than a few compounds (scope tables), give a table and
write it with `export` (CSV for spreadsheets, SDF for ChemDraw/registration
systems).

## Pitfalls

- Stereo: wedges and hashes matter; if a figure shows none, do not invent
  stereo. Report unassigned centres instead.
- Aromatic heterocycles: put the H on the right nitrogen ([nH]); check
  tautomers drawn explicitly.
- Ring fusions and bridged/cage systems are the most common misreads: count
  atoms in every ring of the redraw against the original.
- Counter-ions and salts are separate fragments (`.`): keep them if drawn.
- Superscripts/subscripts and atom labels in low-resolution scans are easy to
  misread (Cl vs C1, O vs 0, NH vs N): zoom further before guessing.
- Polymers, organometallics with multicentre bonds, and Markush structures do
  not fit plain SMILES well: describe them in words with partial SMILES.
