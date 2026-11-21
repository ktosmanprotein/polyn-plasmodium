# polyn-plasmodium

Analysis code for:

> Osman KT. *Asparagine repeats in malaria parasites are the protein face of a genome-wide shift towards tandem repeats.* Submitted to eLife (2026).

The derived data produced by this code (tract coordinates and ages, placements, reconstructions, Pf8 summaries, repeat counts, figure source data) are archived separately on Zenodo: **doi:10.5281/zenodo.DATA_DOI**.

## What the code does

| Step | Scripts | Paper section |
|---|---|---|
| Proteomes, low-complexity regions, masking | `fetch_sequences.py`, `call_tracts.py`, `summarise_lcr.py`, `mask_lcr.py`, `seg_check.py` | Methods: low-complexity regions |
| Orthology and structures | `build_orthogroups.py`, `fetch_structures.py`, `align_orthogroups.py`, `og_tracts.py` | Methods: orthology and structures |
| Structure-anchored placement and insertion sites | `analyse_placement.py`, `insertion_site_null.py`, `insertion_site_dssp.py`, `dssp_lite.py`, `seq_anchor_placement.py` | Results: Figure 2 |
| Segment length, seeding, N versus K | `segment_length_check.py`, `compositional_tests.py` | Results: Figure 3 |
| Dating, ancestral reconstruction, codon decay | `date_tracts.py`, `ancestral_reconstruction.py`, `codon_decay.py`, `codon_decay_within.py` | Results: Figure 4 |
| Pf8 population analyses | `pf8_zarr_extract.js`, `pop_length_variation.py` | Results: Figure 5 |
| Essentiality and expression | `constraint_vs_age.py` | Results: Figure 6 |
| Genome-wide repeats, repair genes | `genome_repeats.js`, `mmr_check.py` | Results: Figure 1, Figure 6—figure supplement 1 |
| Figures, supplements, source data | `make_figures_v2.py`, `make_supplements_sourcedata.py` | All figures |

Shared helpers are in `polylib.py` and `structlib.py`.

## Inputs (public, not redistributed)

- PlasmoDB / ToxoDB release 71: proteins, CDS, genomes and GFF for 13 *Plasmodium* species and *T. gondii* ME49.
- AlphaFold DB version 6 models.
- MalariaGEN Pf8 (`Pf8.zarr.zip`). The JavaScript extractors read it by HTTP range requests from a browser console on the Pf8 release site; see the header of `pf8_zarr_extract.js`.

## Requirements

Python ≥ 3.9 with NumPy, SciPy, Matplotlib and Pillow (`pip install -r requirements.txt`); MMseqs2 (≥ 15); Foldseek; fLPS 2. `dssp_lite.py` is a NumPy re-implementation of DSSP and needs nothing else.

## Citation

Please cite the paper and the archived code release (see `CITATION.cff`).

## Licence

MIT (see `LICENSE`).
