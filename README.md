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
> contain real, tested logic; the external-tool rules carry runnable command
> templates. The Python screens run today (see
> [Running the screens standalone](#running-the-screens-standalone)); the full
> Snakemake DAG runs once the external tools (MAFFT, IQ-TREE, HMMER, TreeShrink,
> GeneRax, PhylteR …) and your marker/reference resources are installed — see
> [Installation](#installation).

---

## Contents

- [Quickstart](#quickstart)
- [How it works (and why)](#how-it-works-and-why)
- [Installation](#installation)
- [Preparing your inputs](#preparing-your-inputs)
- [Configuring a run](#configuring-a-run)
- [Running the pipeline](#running-the-pipeline)
- [Understanding the outputs](#understanding-the-outputs)
- [Running the screens standalone](#running-the-screens-standalone)
- [Pipeline stages and scripts](#pipeline-stages-and-scripts)
- [Extending MarkerVet](#extending-markervet)
- [Troubleshooting / FAQ](#troubleshooting--faq)
- [References](#references)
- [License](#license)

---

## Quickstart

```bash
# 1. clone
git clone https://github.com/Tancata/markervet.git
cd markervet

# 2. create the tool environment (see Installation for alternatives)
mamba env create -f envs/markervet.yaml
conda activate markervet

# 3. drop your per-marker FASTAs + reference tree into place
#    (edit the paths in config/prokaryote_deep.yaml to match)
mkdir -p input/prokaryote_deep/markers
cp /path/to/*.faa input/prokaryote_deep/markers/
cp /path/to/gtdb_reference.tree resources/gtdb/gtdb_reference.tree

# 4. dry-run to inspect the plan, then run
snakemake -n  --configfile config/prokaryote_deep.yaml
snakemake -c8 --configfile config/prokaryote_deep.yaml

# 5. read the verdict
column -t -s$'\t' results/prokaryote_deep/aggregate/scorecard.tsv
open  results/prokaryote_deep/report/markervet_report.html
```

The markers that survive are in `results/prokaryote_deep/outputs/` as cleaned
alignments, collapsed gene trees, and a concatenated supermatrix + partition
file.

---

## How it works (and why)

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
keeping it. Loosen the policy by raising a threshold, or by setting a
criterion's `enabled: false` to remove that vote entirely.

---

## Installation

MarkerVet needs a Python/R stack (for the screens and the workflow engine) plus
several external phylogenetics tools (for alignment, trimming, tree inference,
HMM search, and reconciliation).

### Option A — one conda/mamba environment (recommended)

```bash
mamba env create -f envs/markervet.yaml   # or: conda env create -f ...
conda activate markervet
```

This installs everything: `snakemake`, `mafft`, `clipkit`, `bmge`, `iqtree`
(provides `iqtree2`), `hmmer`, `treeshrink`, `generax`, `ale`, the Python stack
(`numpy`, `scipy`, `pyyaml`, `dendropy`, `ete3`, `pandas`), and the R stack
(`r-base`, `r-phylter`, `r-ape`, `r-optparse`).

### Option B — Snakemake-managed per-rule envs

Keep only Snakemake in your base environment and let it build the tool
environment on first run:

```bash
mamba create -n smk -c conda-forge -c bioconda snakemake-minimal
conda activate smk
snakemake -c8 --use-conda --conda-frontend mamba \
          --configfile config/prokaryote_deep.yaml
```

### Verifying the install

```bash
# tools on PATH
for t in mafft clipkit iqtree2 hmmsearch run_treeshrink.py generax; do
    command -v "$t" >/dev/null && echo "ok   $t" || echo "MISSING $t"
done
# python + R libs
python -c "import numpy, scipy, yaml, dendropy, ete3; print('python deps ok')"
Rscript -e 'library(phylter); library(ape); library(optparse); cat("R deps ok\n")'
```

Pins in `envs/markervet.yaml` are intentionally loose (`>=`) so the scaffold
resolves on current bioconda. **Tighten them to exact versions for a
reproducible production run.**

---

## Preparing your inputs

All inputs are referenced by path from the regime YAML. The default layout the
presets expect:

```
markervet/
├── input/
│   └── prokaryote_deep/
│       ├── markers/           # {marker}.faa  (extract: false, the default)
│       └── proteomes/         # {taxon}.faa   (only if extract: true)
├── resources/
│   └── gtdb/
│       ├── gtdb_reference.tree # reference topology, pruned to your sampling
│       └── hmms/               # {marker}.hmm (only if extract: true)
└── config/
    ├── prokaryote_deep.yaml
    └── prokaryote_deep.taxa.tsv
```

You can put files anywhere — just point the config at them.

### 1. Markers — two ways to supply them

**(a) Supplied per-marker FASTAs** (`markers.extract: false`, the default).
One **unaligned** protein FASTA per marker, already reduced to your taxon
sampling, named `{marker}{fasta_ext}` (e.g. `ar53_arCOG00780.faa`). The marker
name is the filename without the extension and becomes the wildcard everywhere.

**Critical:** the FASTA header (first whitespace-delimited token) *is* the taxon
name, and it must be **identical across every marker and in the reference
tree** — that is how MarkerVet lines taxa up across genes.

```text
>GB_GCA_000010565.1
MSKVLIVGAGPAGL...
>RS_GCF_000008665.1
MSKILAVGAGPAGL...
```

**(b) De novo extraction** (`markers.extract: true`). Supply per-taxon
proteomes as `{taxon}.faa` in `markers.proteomes_dir`, and one profile HMM per
marker as `{marker}.hmm` in `markers.hmm_dir`. MarkerVet runs `hmmsearch`, then
enforces single-copy: it keeps the best-scoring hit per taxon and **drops a
taxon whose second-best hit scores within `paralog_margin` bits of the best**
(ambiguous paralogy). Sequences are relabelled to the taxon name on the way out,
so downstream leaf names are uniform. Header convention for extraction is
`taxon|...`; adjust `taxon_of()` in `rules/extraction.smk` if yours differs.

### 2. Reference tree (required)

A trusted topology in Newick, `reference_tree`, pruned to (a superset of) your
taxon sampling — e.g. the GTDB reference tree, or a curated three-domain tree.
Used by the congruence stage (RF + quartet + gCF/sCF). Only the shared taxa are
compared, so extra tips in the reference are fine.

### 3. Taxon table (optional)

A TSV, `taxon_table`, whose **first column is the leaf name** and whose other
columns are arbitrary taxonomy fields. Needed **only** if a clade definition
targets a `column:` instead of matching the leaf name directly. Examples ship as
`config/*.taxa.tsv`.

```text
taxon	domain	phylum	class
GB_GCA_000010565.1	d__Archaea	p__Nanoarchaeota	c__Nanoarchaeia
RS_GCF_000008665.1	d__Bacteria	p__Cyanobacteria	c__Cyanophyceae
```

---

## Configuring a run

Everything is driven by the regime YAML you pass with `--configfile`. Two
presets ship in `config/`:

| | `prokaryote_deep.yaml` | `three_domain.yaml` |
|---|---|---|
| Markers | GTDB `ar53` + `bac120` | pan-domain universal subset |
| Reference | GTDB reference tree | three-domain reference |
| Recoding track | off | **on** (Dayhoff-6) |
| Compositional / saturation cuts | standard | **stricter** |
| Support gate | strict (few, dense markers) | **looser** (few, sparse markers) |
| Model | LG+C20+G (PMSF) | LG+C60+G (PMSF) |
| Monophyly clades | d__Archaea, d__Bacteria, DPANN, Asgard, Patescibacteria (CPR), example phylum/class | Archaea, Bacteria, Eukaryota, Asgard, TACK, DPANN |

Copy a preset to start your own: `cp config/prokaryote_deep.yaml config/mystudy.yaml`.

### Full configuration reference

| Key | Meaning |
|---|---|
| `regime`, `workdir` | run label; output root (keep distinct per regime so runs don't clobber) |
| `markers.extract` | `false` = use supplied FASTAs; `true` = hmmsearch extraction |
| `markers.fasta_dir`, `markers.fasta_ext` | where per-marker FASTAs live / their extension |
| `markers.proteomes_dir`, `markers.hmm_dir` | proteomes and HMMs (extraction mode) |
| `markers.hmm_list` | explicit marker list (path or inline YAML list); `null` = discover all |
| `markers.hmmsearch_evalue`, `markers.paralog_margin` | hmmsearch E-value; bits within which a paralog makes a taxon ambiguous |
| `taxon_table` | optional per-taxon annotation TSV |
| `reference_tree` | trusted reference topology (Newick) |
| `align.tool` / `align.mode` | `mafft`; `auto` or `linsi` (accurate L-INS-i) |
| `trim.tool` | `clipkit` (with `clipkit_mode`, e.g. `smart-gap`) or `bmge` (with `bmge_entropy`) |
| `artefact.compositional.enabled` / `.alpha` | chi-square composition test; fail if *p* < α |
| `artefact.saturation.enabled` / `.min_slope` | p-dist~patristic slope; fail below |
| `artefact.treeshrink.enabled` / `.quantile` | rogue-tip removal aggressiveness |
| `gene_trees.model` | site-heterogeneous model, e.g. `LG+C20+G` / `LG+C60+G` |
| `gene_trees.pmsf` | two-step PMSF approximation of the mixture (fast) |
| `gene_trees.ufboot`, `.alrt`, `.threads` | UFBoot replicates; SH-aLRT replicates; threads |
| `collapse.ufboot_min`, `.alrt_min` | support thresholds for branch collapse |
| `collapse.logic` | `either` (collapse if it fails either gate — strict) or `both` (lenient) |
| `congruence.enabled` / `.support_gate` | run congruence; support at/above which a branch counts as conflict |
| `congruence.rf_max`, `.quartet_max` | normalised RF / quartet ceilings vs reference |
| `congruence.gcf_min` | optional gene-concordance-factor floor (`null` = ignore) |
| `monophyly.support_gate` | collapse branches below this before scoring clades |
| `monophyly.recovery_min` | drop if `clade_recovery_score` falls below this |
| `monophyly.clades[]` | list of `{name, pattern, column?}` clade definitions |
| `outliers_hgt.phylter.enabled` / `.taxon_cell_frac_max` | PhylteR vote; max fraction of taxon cells it may remove before drop |
| `outliers_hgt.reconciliation.enabled` / `.tool` / `.transfers_max` | reconciliation vote; `generax` or `ale`; transfer budget |
| `recoding.enabled` / `.scheme` | AA recoding track; `dayhoff6` \| `sr4` \| `sr6` |
| `policy.mode` | `aggressive` (drop on any failure) |

**Clade definitions** are the heart of the monophyly stage. Each entry:

```yaml
monophyly:
  clades:
    - {name: d__Archaea,   pattern: "d__Archaea"}        # regex on leaf name
    - {name: Eukaryota,    pattern: "Eukaryota|euk__"}
    - {name: c__Alphaproteobacteria, pattern: "c__Alphaproteobacteria", column: class}
```

`pattern` is a Python regex matched against the leaf name, or — if `column` is
given — against that column of `taxon_table`. First match in list order wins, so
list more specific clades first if they can overlap.

---

## Running the pipeline

All commands are run from the repo root with the environment active.

```bash
# Dry-run: print the plan without executing (always do this first).
snakemake -n --configfile config/prokaryote_deep.yaml

# Visualise the DAG.
snakemake --dag --configfile config/prokaryote_deep.yaml | dot -Tsvg > dag.svg

# Full run on 8 cores.
snakemake -c8 --configfile config/prokaryote_deep.yaml

# Let Snakemake manage the conda env per rule.
snakemake -c8 --use-conda --conda-frontend mamba \
          --configfile config/three_domain.yaml
```

### Useful targets and flags

```bash
# Stop after the scorecard (skip building supermatrix/report).
snakemake -c8 --configfile config/prokaryote_deep.yaml \
          results/prokaryote_deep/aggregate/scorecard.tsv

# Screen a single marker end-to-end (any per-marker output works as a target).
snakemake -c4 --configfile config/prokaryote_deep.yaml \
          results/prokaryote_deep/monophyly/ar53_arCOG00780.summary.tsv

# Re-run everything downstream of a config change.
snakemake -c8 --configfile config/prokaryote_deep.yaml --forceall

# Keep going after a single marker fails, and retry flaky tool calls.
snakemake -c8 --configfile config/prokaryote_deep.yaml --keep-going --retries 2
```

### Cluster / HPC execution

Use a Snakemake **profile** or the generic cluster submission flag; every rule
declares `threads`, so per-job resources map cleanly onto a scheduler.

```bash
# With a profile (recommended): create ~/.config/snakemake/slurm/config.yaml, then
snakemake --profile slurm --configfile config/prokaryote_deep.yaml

# Or the classic form:
snakemake -j100 --configfile config/prokaryote_deep.yaml \
          --cluster "sbatch -c {threads} --mem=8G -t 04:00:00"
```

The pipeline fans out per marker, so it parallelises trivially across a cluster;
the only whole-set step is PhylteR (stage 7), which needs every gene tree at once.

---

## Understanding the outputs

Everything lands under the regime's `workdir` (e.g. `results/prokaryote_deep/`).
Intermediate directories (`align/`, `trim/`, `artefact/`, `genetrees/`,
`congruence/`, `monophyly/`, `outliers_hgt/`) hold the per-stage evidence; the
deliverables are:

```
results/<regime>/
├── aggregate/
│   ├── scorecard.tsv            # ← the verdict: keep/drop + reasons + stats
│   ├── effective_config.yaml    # the fully-resolved config this run used
│   └── markers.txt              # the marker list screened
├── outputs/
│   ├── alignments/{marker}.faa  # cleaned, trimmed alignments (passers only)
│   ├── gene_trees/{marker}.nwk  # collapsed passer gene trees → coalescence
│   ├── supermatrix.fasta        # concatenated passers → concatenation
│   └── supermatrix.partitions   # per-marker charsets for the supermatrix
└── report/
    ├── markervet_report.html    # keep/drop summary + drop-reason tally +
    └── markervet_report.md      #   clade × marker recovery matrix
```

### The scorecard

`scorecard.tsv` is one row per marker. The first columns are the decision; the
rest are the key statistic behind each criterion (handy for tuning thresholds).

| Column | Meaning |
|---|---|
| `marker` | marker name |
| `decision` | `KEEP` or `DROP` |
| `n_failed` | number of failed criteria (0 ⇒ KEEP) |
| `reasons` | `;`-joined failure reasons with the offending value, e.g. `saturated(slope=0.2);excess_transfers(10>6)` |
| `comp_p` | chi-square composition *p*-value |
| `sat_slope` | saturation regression slope |
| `rf_norm`, `quartet_norm` | normalised RF / quartet distance vs reference |
| `clade_recovery_score`, `intruders`, `outliers` | monophyly headline + counts |
| `phylter_frac_removed` | fraction of this marker's taxon cells PhylteR removed |
| `transfers` | reconciliation-inferred transfer count |

Reasons are fail-closed: a `missing:<criterion>` reason means that screen's
result file was absent, and the marker was dropped rather than trusted.

```bash
# quick tallies
column -t -s$'\t' results/<regime>/aggregate/scorecard.tsv     # pretty-print
cut -f2 results/<regime>/aggregate/scorecard.tsv | sort | uniq -c   # KEEP/DROP counts
awk -F'\t' '$2=="KEEP"{print $1}' results/<regime>/aggregate/scorecard.tsv  # passer names
```

### The report

`markervet_report.html` (and `.md`) summarises: the KEEP/DROP split, a tally of
*why* markers were dropped (which lets you see whether one failure mode
dominates), the full decision table, and the **clade × marker recovery matrix**
— a scientific result in its own right, showing which markers recover which deep
clades.

---

## Running the screens standalone

The pure-Python/R screens are ordinary CLI scripts — useful for debugging,
teaching, or applying one test outside the workflow. Each prints its verdict to
stderr and writes a one-row TSV (`--out-tsv`). Examples (the first three need
only `numpy`/`scipy`/`dendropy`; the monophyly/collapse scripts need `ete3`):

```bash
# Compositional homogeneity (chi-square) on one alignment.
python scripts/compositional_test.py \
    --alignment aln.faa --marker myMarker --alpha 0.05 --out-tsv comp.tsv

# Saturation slope (needs an alignment + a tree with branch lengths).
python scripts/saturation.py \
    --alignment aln.faa --tree gene.treefile --marker myMarker \
    --min-slope 0.55 --out-tsv sat.tsv

# Normalised RF + quartet distance vs a reference.
python scripts/tree_distances.py \
    --gene-tree gene.nwk --reference ref.nwk --marker myMarker \
    --rf-max 0.35 --quartet-max 0.30 --out-tsv cong.tsv

# Collapse branches failing the IQ-TREE aLRT/UFBoot gate.
python scripts/collapse_support.py \
    --tree gene.treefile --out gene.collapsed.nwk --marker myMarker \
    --ufboot-min 95 --alrt-min 80 --logic either

# MonoPhy-style clade recovery (clade defs as a JSON list of {name, pattern}).
python scripts/monophyly_recovery.py \
    --tree gene.collapsed.nwk --marker myMarker --clades-json clades.json \
    --support-gate 80 --recovery-min 0.80 \
    --out-matrix mono.matrix.tsv --out-summary mono.summary.tsv

# PhylteR across a directory of collapsed gene trees.
Rscript scripts/phylter_run.R --tree-dir genetrees/ --tree-ext .collapsed.nwk \
    --out-cells phylter.cells.tsv --out-summary phylter.summary.tsv
```

---

## Pipeline stages and scripts

| # | Stage | Rule file | What it does |
|---|---|---|---|
| 1 | extraction | `extraction.smk` | hmmsearch vs marker HMMs + single-copy enforcement (skipped if per-marker FASTAs are supplied) |
| 2 | align_trim | `align_trim.smk` | MAFFT (auto/L-INS-i) → ClipKIT (smart-gap) or BMGE |
| 3 | artefact_screen | `artefact_screen.smk` | compositional chi-square, saturation slope, TreeShrink rogue-tip removal, optional AA recoding |
| 4 | gene_trees | `gene_trees.smk` | IQ-TREE site-heterogeneous model (LG+C20/C60+G, PMSF), UFBoot + aLRT → support-based collapse |
| 5 | congruence | `congruence.smk` | IQ-TREE gCF/sCF + normalised RF and quartet distance vs reference |
| 6 | monophyly | `monophyly.smk` | MonoPhy-style intruder/outlier counts + per-clade recovery, clade × marker matrix |
| 7 | outliers_hgt | `outliers_hgt.smk` | PhylteR (gene × taxon outliers) + GeneRax/ALE reconciliation (excess transfers) |
| 8 | aggregate | `aggregate.smk` | apply the aggressive policy → `scorecard.tsv` (a Snakemake checkpoint) |
| 9 | outputs | `outputs.smk` | passer alignments, per-marker gene trees, concatenated supermatrix + partitions, HTML/markdown report |

Stage flow (the cleaned alignment `clean.faa`, produced by the artefact screen,
is the canonical input everywhere downstream):

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

## Extending MarkerVet

- **New regime.** Copy a preset (`cp config/prokaryote_deep.yaml
  config/mystudy.yaml`), edit paths, thresholds, and clades, then run with
  `--configfile config/mystudy.yaml`. No code changes needed.
- **New clade to test.** Add a `{name, pattern}` entry under
  `monophyly.clades`; it appears automatically in the recovery matrix and
  becomes part of the drop criterion.
- **Swap a tool.** The external-tool rules are thin shell templates — e.g. point
  `trim.tool` at `bmge`, switch `outliers_hgt.reconciliation.tool` to `ale`, or
  edit the `iqtree2` command in `rules/gene_trees.smk`.
- **Make it more/less aggressive.** Raise a threshold (e.g.
  `artefact.saturation.min_slope`), lower a ceiling (e.g. `congruence.rf_max`),
  or set a criterion `enabled: false` to drop that vote.

---

## Troubleshooting / FAQ

**`snakemake -n` reports "no markers discovered".** The `markers.*` paths in
your configfile don't point at any `*.faa` (or `*.hmm` in extraction mode).
Check `markers.fasta_dir` / `markers.fasta_ext` and that files are actually there.

**Everything gets dropped.** The aggressive policy drops on *any* failure, and
is fail-closed. Read the `reasons` column: a wall of `missing:*` means upstream
stages didn't produce results (usually a tool not on PATH — re-run the install
check). Genuine failures concentrated in one reason mean that threshold is too
strict for your data; loosen it, or set that criterion `enabled: false`.

**Congruence is all `NA` / everything fails congruence.** Fewer than four taxa
are shared between a gene tree and the reference. Confirm leaf names are
*identical* across markers and the reference tree (a trailing `|ncbi_id`, a `.1`
version suffix, or `_` vs space differences will silently break the join).

**A marker's gene tree is empty or IQ-TREE errors.** Very small or gap-only
alignments after trimming. Inspect `results/<regime>/trim/{marker}.trimmed.faa`;
consider a gentler trim (`clipkit_mode`, or higher `bmge_entropy`).

**Can I resume after a crash?** Yes — Snakemake is idempotent. Re-run the same
command; completed jobs are skipped. Use `--rerun-incomplete` if a job was
killed mid-write.

**Where do I change the amino-acid recoding scheme?** `recoding.scheme`
(`dayhoff6` / `sr4` / `sr6`); the class maps live in `rules/artefact_screen.smk`.

**The Python scripts need Python 3.** They target Python ≥3.9 (the repo's older
scripts are Python 2; MarkerVet's are not).

---

## References

MarkerVet orchestrates established methods; please cite the underlying tools you
use:

- **MAFFT** — Katoh & Standley 2013, *MBE*.
- **ClipKIT** — Steenwyk et al. 2020, *PLoS Biol*. **BMGE** — Criscuolo & Gribaldo 2010, *BMC Evol Biol*.
- **TreeShrink** — Mai & Mirarab 2018, *BMC Genomics*.
- **IQ-TREE 2** — Minh et al. 2020, *MBE*; **UFBoot2** — Hoang et al. 2018; **PMSF** — Wang et al. 2018; **concordance factors** — Minh et al. 2020.
- **PhylteR** — Comte et al. 2023, *MBE* (successor to Phylo-MCOA, de Vienne et al. 2012).
- **GeneRax** — Morel et al. 2020, *MBE*; **ALE** — Szöllősi et al. 2013.
- **MonoPhy** (concept) — Schwery & O'Meara 2016, *PeerJ*.
- **GTDB** markers/reference — Parks et al. 2018/2022.
- **Snakemake** — Mölder et al. 2021, *F1000Research*.

---

## License

MarkerVet is released under the **GNU General Public License v3.0** — see
[`LICENSE`](LICENSE). The external tools it orchestrates carry their own
licenses; consult each project when redistributing.

---

*MarkerVet is part of the [tancata/phylo](https://github.com/tancata/phylo)
family of deep-phylogeny tooling.*
