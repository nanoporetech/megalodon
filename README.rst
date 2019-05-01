.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

Megalodon provides "basecalling augmentation", including direct, reference-guided SNP and modified base calling.

Megalodon anchors the information rich neural network to a reference genome. Variants, either modified bases or alternative bases, are then proposed and scored in order to produce highly-accurate reference anchored calls.

Features
--------

- Single command interface

  - Inputs

    - Raw reads
    - Sequence reference

    - Optionally

      - Variants VCF (required for SNP calling)

  - Outputs

    - Basecalls

      - Including basecalls with annotated modified bases
      - Or separate HDF5 modified base scores

    - Mappings (via minimap2 python interface mappy)

      - In either SAM, BAM or CRAM format

    - Modified base calls

      - Per-read modified base calls (either database or text)
      - Aggregated calls (modVCF TODO: add description)

    - Sequence variant calls

      - Per-read sequence variant calls (either database or text)
      - Aggregated calls (VCF)
      - Future additions:

          - Phased read calls
          - Improved phase-aware per-read sequence variants calls
          - Phased VCF output
