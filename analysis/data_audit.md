# Data audit

Checked against the sources, 19 August 2026.

## MegaScale (Tsuboyama 2023)

Zenodo 7992926. Taken: `K50_dG_tables.zip` (450 MB), `AlphaFold_model_PDBs.zip`
(14 MB), `Raw_NGS_count_tables.zip` (243 MB). The HuggingFace mirror is prepared
for ThermoMPNN training and duplicates the Zenodo tables, so it was skipped. The
raw NGS counts were retrieved and are not used by any analysis.

1,868,872 rows across four libraries: Lib1 262,953, Lib2 553,699, Lib3 773,374,
Lib4 278,846. Every row has a parseable `deltaG` and a `dna_seq` in column two,
and every one is in frame. Construct lengths are near enough fixed per library:
132 or 150 nt (Lib1), 216 nt (Lib2, Lib3), 240 nt (Lib4), with 59 rows at 129.

### The DNA is library construction, not a deployment sequence

Codon usage over the 772 recovered wild-type constructs is restricted to two
codons per degenerate residue:

| Residue | Codons observed |
|---|---|
| Leu | CTG 0.83, CTC 0.17 |
| Arg | CGT 0.91, CGC 0.09 |
| Ser | TCT 0.71, TCC 0.29 |
| Ala | GCG 0.71, GCT 0.29 |
| Val | GTT 0.85, GTG 0.15 |
| Gly | GGT 0.74, GGC 0.26 |
| Pro | CCG 0.79, CCA 0.21 |

Leucine has six codons and two appear. This reflects oligo-pool synthesis
constraints rather than a codon optimiser. Designs are padded to constant length
with GS and SAG linkers, 0 to 4 residues at the N terminus and up to 11 at the
C terminus, leaving cores of 39 to 75 residues.

Consequence: the descriptive half of H1, that designs as ordered carry a
measurable liability load, cannot be tested on these sequences. The method
comparison is unaffected, being CDS against CDS by construction.

### Sizing

H4 needs a backbone, so only parents with a model in `AlphaFold_model_PDBs`
count.

| | parents | with singles | single measurements | with doubles | double measurements |
|---|---|---|---|---|---|
| with a backbone | 767 | 730 | 750,482 | 471 | 387,989 |
| de novo designed | 352 | 317 | 306,922 | 171 | 20,449 |
| natural | 415 | 413 | 443,560 | 300 | 367,540 |

Double-mutant coverage on designs is uneven: of 171 designed parents with
doubles, 134 have fewer than ten, 35 have 100 or more and 8 have over 500
(`HHH_rd2_0165` 1,472, `EHEE_rd1_0284` 1,459). Epistasis is workable on those 35.
Wider claims rest on the natural domains, where coverage is far deeper, and that
transfer must be stated.

Composition-matched scrambled variants (`scramble_50%`) are present and cover
part of the D9 control set.

## Rocklin 2017

Science supplementary returns 403 to automated requests; no workaround
attempted. PMC carries the article but not the data tables. The Rocklin lab
GitHub holds analysis pipelines only.

Largely superseded: MegaScale contains 83,728 parents matching the Rocklin
topology grammar (HEEH 75,056, EEHEE 4,141, EHEE 2,601, HHH 1,930, rounds 1 to
6), 252 with a backbone, measured as absolute folding free energy rather than a
proteolysis-derived stability score. Rocklin 2017 uniquely holds the
yeast-display expression readout and the oligo-pool designs, which need a manual
download.

## Garcia, Dixit and Rocklin 2026

*Protein Science* 35(2):e70453, published 20 January 2026, doi:10.1002/pro.70453.
Preprint bioRxiv 2025.07.29.667290 (CC-BY), v1 posted 1 August 2025 and v2 posted
9 August 2025; the copy held in `external/papers` is v1. The preprint carries no
data availability statement and links no supplementary files from the article
page, and the journal version is behind a publisher block, so the 614-design
benchmark was downloaded by hand from the preprint supplement on 20 August 2026.
Cite the journal version and record the preprint as the route to the data.

## PDB de novo subset

2,056 entries, 2,423 polymer entities, median length 100, 571 in the 40 to 120
range. Release dates: 1,334 entities in 2021 or earlier, 1,089 in 2022 or later.
This supplies the H5 contamination stratification.

## Backbones

862 AlphaFold models extracted, 26 to 74 residues, median 47, of which 447 are
designs. These are the working set.

## Consequence for the design

DNA exists but records library synthesis, so the descriptive gene-layer claim in
H1 is dropped and H1 runs on protein-layer features against the control set. The
method comparison is unaffected. H5 is scoped to the protein layer from the
outset, with the release-date stratification available.

The control set as built is four arms: scrambled, composition-matched random, the
excised natural domains from this release, and a whole-protein arm assembled
separately from UniProt. The last two are never pooled. Gene-layer features are
still computed for every arm, but on a most-adapted back-translation rather than
on any deposited DNA, which is why CAI is exactly 1.000 in every row of the
atlas.
