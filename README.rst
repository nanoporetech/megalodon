.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

Megalodon provides "basecalling augmentation", including direct, reference-guided SNP and modified base calling.

Megalodon anchors the information rich neural network basecalling output to a reference genome. Variants, either modified bases or alternative bases, are then proposed and scored in order to produce highly-accurate reference anchored calls.


Installation
------------

Requires taiyaki installation for basecalling backend.

::

    git clone https://git.oxfordnanolabs.local/algorithm/megalodon
    cd megalodon
    pip install .

Getting Started
---------------

Megalodon is accessed via the command line interface, ``megalodon`` command.

::

    # megalodon help (common args)
    megalodon -h
    # megalodon help (all args)
    megalodon --help-long
    
    # example command calling variants and CpG methylation
    megalodon reads/ --taiyaki-model-filename taiyaki/models/mGru_flipflop_remapping_model_r9_DNA.checkpoint \
        --outputs basecalls mappings snps mods --reference reference.fa \
        --snp-filename variants.vcf --snp-calibration-filename lhr_calibration.npz \
        --mod-motif Z CG 0 --devices 0 --processes 8 --verbose-read-progress 3

This command produces the ``megalodon_results`` output directory eith the following files present: ``basecalls.fasta``, ``mappings.bam``, ``per_read_snp_calls.db``, ``snps.vcf``, ``per_read_modified_base_calls.db`` and ``mods.mvcf``.

Simple Command Interface
------------------------

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
