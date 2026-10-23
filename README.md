# Ubiquitin contact patterns on a disordered polystyrene nanoparticle

Reproducibility materials for **Recurring Ubiquitin Contact Patterns on a Disordered Polystyrene Nanoparticle**, by **Benjamin Borow and Ashley Z. Guo** (manuscript intended for *Langmuir*).

This package contains the numerical tables and simulation-derived molecular images used to present the first interaction/mechanism study. Its fixed scientific dataset is `AIM1-MECHANISM-20261016-v1`, with descriptive surface addendum `AIM1-SURFACE-20261021-v1`. It contains no results from the separate affinity campaign or enhanced-sampling benchmark.

## What the evidence supports

An irregular polystyrene particle accommodates several ubiquitin contact patterns. The sampled interfaces share recurring chemical, direct-attraction and hydration features while retaining much of the sampled protein core. The three native preparation histories are correlated, selected histories; the orientation-release families are a different classification. Neither is an equilibrium population distribution. The B distance profile is conditional on its sampled preparation and is not a preparation-independent adsorption affinity. Direct interaction energies and regional hydration terms are not an additive decomposition of binding free energy. Surface-map differences persist over the sampled interval but do not demonstrate permanent, causally preparation-induced particle deformation.

## Reproduce figures from the supplied tables

Eight larger figure and coordinate files are distributed as exact binary parts. The first command restores their original paths, checks every part and final file hash, and refuses to overwrite changed files. No scientific bytes are altered. Keep `large_file_parts/` and `LARGE_FILE_TRANSPORT.json` together with the repository.

Use Python 3.12 (the verified version) and a virtual environment:

```bash
python restore_large_files.py
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python verify_package.py
python reproduce_figures.py
python make_toc.py
python figure_revision/reproduce_figures.py
```

The first plotting command rebuilds the six main figures and the conditional PMF/interaction-energy comparison (supporting Figure S11) from the supplied numerical tables. The supplementary command rebuilds Figures S10 and S12 and independently checks the clustering and ring-geometry summaries from supplied mean occupancy maps and tabulated measurements. The manuscript-ready PDFs are retained in `manuscript_figures/`. Other supporting figures are supplied as fixed assets with their available source tables; this package does **not** claim that a command reruns every supporting analysis from raw trajectories. See `FIGURE_COVERAGE.md` for the exact coverage.

The main plotting code writes `reproduced_figures/`, `reproduced_tables/`, `supporting_figures/`, `review/` and generated JSON summaries. It does not alter `inputs/`. Molecular Figure 1 uses three archived, simulation-derived renders. Matching coordinates and relative-path ChimeraX scripts are in `molecular_inputs/`; the scripts may be run from that directory with ChimeraX installed separately. Molecular rerendering is not needed for the Python figure build. Small rendering differences can occur across fonts/library versions; numerical tables and semantic checks are the reproducibility target, not PDF creation timestamps.

## Files

| Location | Content |
| --- | --- |
| `inputs/`, `INPUT_MAP.json` | Frozen numerical inputs used directly by the main plotting code. Logical input identifiers retain their original analysis names. |
| `data/plotted/` | Plotted tables, contact/chemistry/structure summaries and quantitative comparisons. |
| `data/supporting/` | Additional numerical supporting data retained in the frozen scientific review package. |
| `molecular_inputs/` | Three medoid coordinate snapshots, figure renders, relative-path ChimeraX commands and residue mapping. |
| `figure_revision/` | Repaired S10/S12 plotting code, original input tables/maps and independently recomputed clustering/ring results. |
| `surface/` | Surface comparisons and sensitivity tables and three mean occupancy maps. |
| `manuscript_figures/` | Figure PDFs synchronized with the revised manuscript. |
| `SOURCE_PROVENANCE.json` | Original and public-file fingerprints; any locator-only sanitization is disclosed. |
| `PACKAGE_MANIFEST.json` | SHA-256 hashes of all distributed files. |
| `LARGE_FILE_TRANSPORT.json`, `large_file_parts/` | Hash-pinned binary parts and exact original-file restoration map. |

Full trajectories, restart files, MD engine/force-field distributions and raw GIST grids are **not included**. This repository supports reproduction and inspection of the reported figure-level results, not regeneration of every trajectory-level measurement. Raw data remain in the authors' preservation archives; no public raw-data URL is asserted here. The raw-data access/deposition arrangement must be described accurately in the submitted paper.

## Citation, rights and release state

Use the authors/title in `CITATION.cff`; no journal DOI is assigned in this release. External articles and third-party software are referenced, not redistributed. See `RIGHTS.md`. Release **2026.10.23** is publicly available at [https://github.com/BBoiSizzle/ubiquitin-polystyrene-mechanisms](https://github.com/BBoiSizzle/ubiquitin-polystyrene-mechanisms). No archival DOI is assigned. A GitHub repository by itself does not supply a persistent archival DOI.
