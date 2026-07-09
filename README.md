# MarkerVet

A config-driven **Snakemake** pipeline that screens phylogenetic marker genes
for **congruence** and **horizontal gene transfer (HGT)** so that only
phylogenetically well-behaved markers reach concatenation and coalescence
analyses of **deep** relationships:

* **(a) prokaryote domain-level trees** — Archaea + Bacteria, and
* **(b) the placement of eukaryotes** relative to Archaea and Bacteria.

Given candidate markers sampled across a fixed taxon set, MarkerVet returns, per
marker, a **keep/drop decision** plus **cleaned, trimmed alignments** and
**collapsed gene trees** ready for downstream supermatrix and gene-tree
(coalescence) methods.

> **Status:** this is a clean, fully-wired *scaffold*. The pure-Python/R steps
> contain real logic; the external-tool rules carry runnable command templates.
> It will not run end-to-end until the external tools (MAFFT, IQ-TREE, …) and
> the marker/reference resources are installed — see `envs/markervet.yaml`.

---

## Why these design choices (deep phylogeny is different)

At domain-level depths, single-gene signal is faint and easily corrupted, so the
pipeline is built around four principles:

1. **Collapse low-support branches *before* any tree comparison.** A raw
   single-gene tree at this depth is mostly estimation error; comparing
   uncollapsed trees measures noise, not conflict. Everything downstream of gene
   inference works on support-collapsed trees (`scripts/collapse_support.py`).
2. **Artefact screen *precedes* HGT calling.** Compositional attraction and
   substitution saturation — the classic eukaryote-placement failure modes —
   are screened first (`compositional_test.py`, `saturation.py`, TreeShrink) so
   they are never mislabelled as transfer.
3. **Two independent HGT/outlier votes.** PhylteR (broad, DISTATIS-based
   gene × taxon outlier detection — the modern successor to Phylo-MCOA) *and*
   reconciliation (GeneRax/ALE transfer localisation).
4. **Monophyly recovery is both a filter and a result.** The clade × marker
   recovery matrix is a headline scientific output as well as a drop criterion.

### Drop policy — aggressive by default

For concatenation, even a little incongruent data biases the tree, so MarkerVet
uses **conservative-keep / aggressive-drop**: a marker is **dropped if it fails
ANY enabled criterion** —

* compositional homogeneity (chi-square) fail,
* saturation slope below threshold,
* PhylteR gene-level outlier **or** > *X*% taxon cells removed,
* reconciliation transfers above threshold,
* high-support RF / quartet conflict with the reference,
* clade-recovery score below threshold.

Every threshold is exposed in the regime YAML. The aggregator is also
*fail-closed*: a missing screen result drops the marker rather than silently
keeping it.

---

## Regime presets (`config/`)

| | `prokaryote_deep.yaml` | `three_domain.yaml` |
|---|---|---|
| Markers | GTDB `ar53` + `bac120` | pan-domain universal subset |
| Reference | GTDB reference tree | three-domain reference |
| Recoding track | off | **on** (Dayhoff-6) |
| Compositional / saturation cuts | standard | **stricter** |
| Support gate | strict (few, dense markers) | **looser** (few, sparse markers) |
| Model | LG+C20+G (PMSF) | LG+C60+G (PMSF) |
| Monophyly clades | d__Archaea, d__Bacteria, DPANN, Asgard, Patescibacteria (CPR), example phylum/class | Archaea, Bacteria, Eukaryota, Asgard, TACK, DPANN |

Named clades are defined by a regex matched against the leaf name (or a column
of the optional taxon table). Edit them freely — they are both the filter and
the reported recovery matrix.

---

## Pipeline stages (`rules/*.smk`)

| # | Stage | Rule file | What it does |
|---|---|---|---|
| 1 | extraction | `extraction.smk` | hmmsearch vs marker HMMs + single-copy enforcement (skipped if per-marker FASTAs are supplied) |
| 2 | align_trim | `align_trim.smk` | MAFFT (auto/L-INS-i) → ClipKIT (smart-gap) or BMGE |
| 3 | artefact_screen | `artefact_screen.smk` | compositional chi-square, saturation slope, TreeShrink rogue-tip removal, optional AA recoding |
| 4 | gene_trees | `gene_trees.smk` | IQ-TREE site-heterogeneous model (LG+C20/C60+G, PMSF), UFBoot + aLRT → support-based collapse |
| 5 | congruence | `congruence.smk` | IQ-TREE gCF/sCF + normalised RF and quartet distance vs reference |
| 6 | monophyly | `monophyly.smk` | MonoPhy-style intruder/outlier counts + per-clade recovery, clade × marker matrix |
| 7 | outliers_hgt | `outliers_hgt.smk` | PhylteR (gene × taxon outliers) + GeneRax/ALE reconciliation (excess transfers) |
| 8 | aggregate | `aggregate.smk` | apply the aggressive policy → `scorecard.tsv` |
| 9 | outputs | `outputs.smk` | passer alignments, per-marker gene trees, concatenated supermatrix + partitions, HTML/markdown report |

Stage flow (cleaned alignment `= clean.faa`, produced by the artefact screen, is
the canonical input everywhere downstream):

```
markers/{m}.faa
   └─ MAFFT ─ align/{m}.aln.faa
        └─ ClipKIT/BMGE ─ trim/{m}.trimmed.faa
             ├─ quicktree ─┐
             └─────────────┴─ TreeShrink ─ artefact/clean/{m}.clean.faa
                  ├─ compositional_test.py ─ artefact/composition/{m}.tsv
                  ├─ saturation.py          ─ artefact/saturation/{m}.tsv
                  └─ IQ-TREE (+PMSF) ─ genetrees/{m}.treefile
                       └─ collapse_support.py ─ genetrees/{m}.collapsed.nwk
                            ├─ tree_distances.py     ─ congruence/{m}.tsv
                            ├─ monophyly_recovery.py ─ monophyly/{m}.summary.tsv
                            ├─ phylter_run.R (all m) ─ outliers_hgt/phylter/summary.tsv
                            └─ GeneRax/ALE           ─ outliers_hgt/recon/{m}.tsv
                                 └─ aggregate_scorecard.py ─ aggregate/scorecard.tsv  [checkpoint]
                                      └─ outputs/{supermatrix, gene_trees, report}
```

---

## Scripts (`scripts/`) — the real logic

| Script | Method |
|---|---|
| `compositional_test.py` | Chi-square test of taxon × amino-acid homogeneity; fail if *p* < α. |
| `saturation.py` | Regression slope of observed p-distance on patristic distance (dendropy); fail if slope < threshold. |
| `collapse_support.py` | Parse IQ-TREE `aLRT/UFBoot` internal labels and collapse branches failing either/both gates (ete3). |
| `tree_distances.py` | Normalised Robinson-Foulds + quartet distance vs the reference on the shared taxon set (dendropy). |
| `monophyly_recovery.py` | MonoPhy-style intruder/outlier counts + per-clade recovery fraction, support-gated; emits clade × marker rows and a per-marker `clade_recovery_score`. |
| `phylter_run.R` | PhylteR DISTATIS gene × taxon outlier detection across all markers; per-marker cell-removal summary. |
| `aggregate_scorecard.py` | Applies the aggressive policy with per-reason drop logic → `scorecard.tsv`. |

---

## Inputs

Two ways to supply markers (per regime YAML, `markers:` block):

* **Supplied per-marker FASTAs** (`extract: false`, the default): one unaligned
  protein FASTA per marker at `markers.fasta_dir/{marker}{fasta_ext}`, already
  reduced to your taxon sampling. Extraction is skipped.
* **De novo extraction** (`extract: true`): per-taxon proteomes in
  `markers.proteomes_dir` plus one HMM per marker in `markers.hmm_dir`.
  hmmsearch + single-copy enforcement build the per-marker FASTAs.

Also required: a **reference tree** (`reference_tree`) pruned to the sampling,
and — optionally — a **taxon table** (`taxon_table`) if any clade definition
targets a taxonomy `column` rather than the leaf name. Leaf names must be
identical across every marker FASTA and the reference tree.

---

## Running

```bash
# create the environment
mamba env create -f envs/markervet.yaml
conda activate markervet

# dry-run to inspect the DAG for a regime
snakemake -n --configfile config/prokaryote_deep.yaml

# run (Snakemake-managed conda envs optional)
snakemake --cores 8 --configfile config/prokaryote_deep.yaml
snakemake --cores 8 --configfile config/three_domain.yaml --use-conda
```

Key outputs land under the regime's `workdir` (e.g. `results/prokaryote_deep/`):

* `aggregate/scorecard.tsv` — per-marker decision + reasons + key statistics
* `outputs/alignments/*.faa` — cleaned, trimmed alignments of passers
* `outputs/gene_trees/*.nwk` — collapsed passer gene trees (coalescence input)
* `outputs/supermatrix.fasta` + `outputs/supermatrix.partitions` — concatenation input
* `report/markervet_report.html` (and `.md`) — keep/drop summary, drop-reason
  tally, and the clade × marker recovery matrix

---

## Tuning the aggressiveness

Loosen or tighten the policy entirely from the regime YAML — e.g. raise
`artefact.saturation.min_slope` or lower `congruence.rf_max` to drop more; set a
criterion's `enabled: false` to remove that vote. Because the policy drops on
*any* failure, disabling a screen is the primary lever for a more permissive
run. Consider the three-domain preset's looser support gate and recoding track a
worked example of adapting the defaults to sparse, deep data.
