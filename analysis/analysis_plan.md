# Analysis plan

Fixed before the atlas was computed.

## Feature panel

Motif definitions come from the Eukaryotic Linear Motif resource, classes file
version 1.4 downloaded 19 August 2026, which supplies curated regular
expressions for 353 classes including 33 degron classes. Using ELM rather than
transcribing motifs from the primary papers means the definitions are versioned,
citable and identical to what other groups scan with.

### Protein layer

Degron and motif counts, each computed twice, once as a raw count and once
weighted by the accessibility of the matching residues:

1. C-end degron hits, ELM `DEG_Cend_*` (DCAF12, FEM1A/C, FEM1B, KLHDC2, TRIM7)
2. N-end degron hits, ELM `DEG_Nend_*` (Nbox, UBRbox 1 to 4)
3. Internal degron hits, remaining ELM `DEG_*` classes
4. Protease cleavage motifs, ELM `CLV_*`
5. Targeting motifs, ELM `TRG_*`, which covers NLS and NES
6. Post-translational modification sites, ELM `MOD_*`
7. Total degron hits across all `DEG_*`

Composition and biophysics:

8. N-terminal residue class under the N-end rule
9. C-terminal residue identity
10. relative accessibility of the N-terminal residue
11. relative accessibility of the C-terminal residue
12. mean relative accessibility
13. longest exposed hydrophobic run, defined below
14. exposed hydrophobic area, defined below
15. count of exposed hydrophobic patches
16. GRAVY
17. net charge at pH 7.4
18. isoelectric point
19. free cysteine count
20. low-complexity fraction
21. helix fraction
22. strand fraction
23. coil fraction

### Gene layer

24. initiation-window folding free energy, minus 4 to plus 37 about the start codon
25. whole-transcript folding free energy
26. codon adaptation index
27. mean codon-pair score
28. overall GC fraction
29. maximum GC in a 20 base window
30. maximum GC in a 100 base window
31. GC over the first 60 bases
32. longest A or T homopolymer
33. longest G or C homopolymer
34. longest exact repeat
35. count of repeats of 20 bases or more
36. internal Shine-Dalgarno hits, AGGAGG allowing one mismatch, outside the initiation window
37. restriction-site hits against the host's vendor panel

Eukaryote-only terms (cryptic splice sites, polyadenylation signals, Kozak
compliance) are computed for the HEK293 host only.

## Thresholds

- exposed: relative accessibility above 0.25
- buried: relative accessibility below 0.15
- hydrophobic: Kyte and Doolittle value above 1.8, which selects A, C, I, L, M, F, V
- exposed hydrophobic run: consecutive residues that are both exposed and hydrophobic; the reported feature is the longest such run
- exposed hydrophobic area: summed absolute accessible area over residues that are both exposed and hydrophobic
- exposed hydrophobic patch: a run of two or more such residues
- low complexity: fraction of positions inside a 12-residue window whose Shannon entropy over residue identity is below 1.5 nats
- accessibility weighting of a motif hit: the mean relative accessibility of its matching residues, so a hit buried in the core contributes near zero and a hit on a free terminus contributes near one
- relative accessibility uses Shrake and Rupley with the theoretical maxima of Tien et al. 2013
- secondary structure uses biotite's P-SEA rather than DSSP, which is unavailable in the container; recorded as a substitution

## Statistics

Primary statistic per feature is the normalised position of the design
distribution on the random-to-natural axis,

    (median_design - median_random) / (median_natural - median_random)

with a bias-corrected and accelerated bootstrap confidence interval over 2,000
resamples. Where the natural and random medians differ by less than the pooled
median absolute deviation the ratio is unstable and the feature is reported as
an effect size only, with the normalised position suppressed.

Secondary statistic is Cliff's delta between designs and each control set, with
a 2,000-resample bootstrap interval.

Correction is Benjamini-Hochberg at a false discovery rate of 0.05 across the
full panel, treated as one family.

Stratification by topology family throughout: HHH, EHEE, EEHEE, HEEH,
hallucination, other.

Every degron and aggregation feature is reported twice, accessibility-weighted
and unweighted. The difference between the two is itself a result.

## Control sets

- natural: *E. coli* BL21(DE3) reference proteins matched to each design on
  length within 10 percent; a second stratum additionally matched on amino-acid
  composition by nearest neighbour in composition space
- scrambled: MegaScale `scramble_50%` variants where available, otherwise a
  within-design permutation of the residues, which preserves composition exactly
  and destroys structure
- random: draws from the aggregate residue composition of the design corpus

Matching is done per design, not per corpus, and the matching record is kept.

## H5 stratifications

Both are mandatory and reported separately.

1. Protein-layer score alone against the full score, reporting the number of
   designs with a recoverable coding sequence. Per the data audit this is
   expected to be zero or near zero, and that number is reported rather than
   glossed.
2. PDB release date, splitting entries released in 2021 or earlier from those
   released in 2022 or later. The split is already computed: 1,334 and 1,089
   polymer entities respectively.

An apparent improvement that does not survive either stratification is not
counted as one.

## Gene-layer scope

Per the data audit, gene-layer features are computed on coding sequences
generated by JANUS and by the baselines, never on the MegaScale as-ordered
sequences. Those are oligo-pool constrained, restricted to two codons per
degenerate residue, and cannot carry a descriptive claim about what designers
order. This has to stay visible in the figure caption.

## Seeds

Global seed 0. Bootstrap seed 1. Control-set sampling seed 2. Permutation-null
seed 3. Recorded in `analysis/config.yaml` and passed explicitly to every
script.

## Is the amino-acid axis worth the paper

Taken 20 August 2026, on the evidence of the exchange rate (`aa_axis.md`), the
delta sweep (`delta_sweep.md`) and the three-arm case study (`case_study.md`).

The question was whether the joint lattice earns its title, or whether the
paper is an exactly-solvable multi-objective codon designer with a measured
amino-acid extension shown to be largely unnecessary. Three measurements decide
it, and they agree.

**A codon-only arm returns exactly zero on every protein-level liability.** Not
approximately zero: 0.000e+00, across 3129 runs per term over seven weights and
447 backbones, on the same hill climb with the residue moves withheld. Against
that, the joint arm clears a low-complexity window outright on 256 of the 257
designs that carry one, and a degron outright on 164 of the 197 that carry one,
for a median 0.21 and 0.22 nats respectively, about half a percent of the shell's
entropy budget.

**Under the complete objective the joint lattice beats the fixed-protein codon DP
on 447 of 447 backbones.** Median gain 1.69 in spread-normalised units at the
operating point, tenth percentile 0.41. There is no subset of backbones on which
the fixed-protein arm wins.

**With the protein held fixed the codon layer's contribution to those liabilities
is exactly zero in the case study too**, in both hosts, while the residue layer
clears degrons, low-complexity content and repeats and improves the initiation
window by more than the codon layer managed.

The title stands and the joint claim leads. Two things constrain how it is
stated.

The amino-acid exchange rate becomes the headline figure rather than the 5' one.
The initiation result stays as its own counterweight: on the term with the
strongest mechanistic evidence and the largest raw effect, 87.6% of the cheap
gain is synonymous, so most of that particular argument is an argument for a good
codon optimiser. The joint lattice earns its place on the terms the codon layer
cannot reach, not on the one it can.

And the axis does not always pay. Between 43% and 56% of designs carry no
low-complexity window and no degron at all, so on those it buys nothing on these
terms. The spread on the delta sweep runs a factor of 34 between the tenth and
ninetieth percentiles. The claim is conditional and the paper states the
condition: when a protein-level liability is present and removable, the codon
layer's achievable gain is zero and the joint lattice's is not, at well under a
nat.

One operating-point result falls out of the same sweep. Half a nat of shell
captures 89% of the achievable gain and one nat captures 99.2%; an unbounded
shell adds under one percent. The delta chosen for the atlas before this was
measured is the right one, and the paper can now say why.

## What this repository does not contain

Figures, the code that draws them, and the manuscript are kept locally and never
enter this repository or any archive made from it. The analysis scripts here
compute and write numbers; every plotting step happens outside.

`restriction.py` used to do both and now writes only its JSON output.
