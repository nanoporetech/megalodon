.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

Megalodon provides "basecalling augmentation" for raw nanopore sequencing reads, including direct, reference-guided SNP and modified base calling.

Megalodon anchors the information rich neural network basecalling output to a reference genome.
Variants, modified bases or alternative canonical bases, are then proposed and scored in order to produce highly-accurate reference anchored modified base or SNP calls.

Detailed documentation for all ``megalodon`` arguments and algorithms can be found on the `megalodon documentation page <https://nanoporetech.github.io/megalodon/>`_.

Installation
------------

Requires `taiyaki <https://github.com/nanoporetech/taiyaki>`_ installation for basecalling backend.

::

    git clone https://git.oxfordnanolabs.local/algorithm/megalodon
    cd megalodon
    pip install .

Getting Started
---------------

Megalodon is accessed via the command line interface ``megalodon`` command.

::

    # megalodon help (common args)
    megalodon -h
    # megalodon help (all args)
    megalodon --help-long

    # example command calling variants and CpG methylation (on GPU devices 0 and 1)
    megalodon raw_fast5s/ --taiyaki-model-filename taiyaki/models/mGru_flipflop_remapping_model_r9_DNA.checkpoint \
        --outputs basecalls mappings snps mods \
        --reference reference.fa --snp-filename variants.vcf \
        --mod-motif Z CG 0 --devices 0 1 --processes 8 --verbose-read-progress 3

This command produces the ``megalodon_results`` output directory containing basecall, mapping, SNP and modified base results.
The format for each output is described below.

Inputs
------

- Raw reads

  - Directory containing raw read FAST5 files
  - By default the directory will be searched recursively for read files (ending in ``.fast5``)
- Reference

  - Genome or transcriptome sequence reference file in FASTA format
- Variants VCF (optional, but required for SNP calling)

  - Megalodon requires a set of candidate variants in order to call SNPs. These should be provided in the VCF format.
  - Only small simple indels (default ``5``) are included in testing. Larger indels can be processed using the ``--max-snp-size`` argument.

    - For larger indels the ``--snp-context-bases`` option may need to be increased.
  - Multiple alternative allele SNPs are not currently supported.

Outputs
-------

- Basecalls

  - Format: FASTA

    - FASTQ format output is not currently available
  - Basecalls with annotated modified bases (TODO describe this output behavior)
  - Or separate HDF5 modified base scores
- Mappings

  - Format: SAM, BAM (default), or CRAM
  - A tab-separated mapping summary is produced

    - Columns: ``read_id``, ``percent_identity``, ``num_aligned_bases``, ``num_matched_bases``, ``num_deleted_bases``, ``num_inserted_bases``

      - ``percent_identity`` is defined as ``num_matched_bases`` / ``num_align_bases``
- Modified Base Calls

  - Per-read modified base calls

    - Per-read SQL DB containing scores at each tested reference location

      - Contains a single ``mods`` table indexed by reference position
    - Tab-delimited output can be produced by adding the ``--write-mods-text`` flag

      - Columns: ``read_id``, ``chromosome``, ``strand``, ``position``, ``score``, ``motif``, ``modified_base``

        - Position is 0-based
        - Motif is as described by ``--mod-motif`` argument

          - If ``--mod-motif`` is not provided, all applicable positions for a modification are tested
  - Aggregated calls

    - Aggregated calls are output in a variant of the VCF format, as no current format allows the output of mulitple types of modifications to the same file.

      - This format treats modified bases as a variant. As opposed to SNP calls (as in VCF format) which output the probability of a particular genotype, this format outputs the estimated proportion of reads modified at the specified genomic location.
- SNP Variant Calls

  - Per-read SNP Calls

    - SQL DB containing scores at each tested reference location

      - Contains a single ``snps`` table indexed by reference position
    - Tab-delimited output can be produced by adding the ``--write-snps-text`` flag

      - Columns: ``read_id``, ``chromosome``, ``strand``, ``position``, ``score``, ``ref_seq``, ``alt_seq``, and ``snp_id``

        - Position is 0-based
  - Aggregated calls

    - Format: VCF
    - VCF file contains ``GT``, ``GQ``, and ``PL`` sample fields
    - Default run mode is diploid. To run in haploid mode, set ``--haploid`` flag.
  - Future additions:

    - Phased read calls
    - Improved phase-aware per-read sequence variants calls
    - Phased VCF output

Computing
---------

Megalodon processes reads from a queue using a pool of workers.
The number of workers is set using the ``--processes`` argument.
Each process is linked to a taiyaki basecalling backend.

In order to use GPU resources the ``--devices`` argument can be set.
If ``--devices`` is set, the taiyaki backends will be distribured evenly over the specified ``--devices``.
In order to control the GPU memory usage, the ``--max_concurrent_chunks`` argument allows a user to restrict the maximum number of chunks to process concurrently (per ``--process``).

The ``--chunk_size`` and ``--chunk_overlap`` arguments allow users to specify read chunking, but signal normalization is always carried out over the entire read.

Compatibility
-------------

The model and calibration files included with megalodon are applicable only to MinION or GridION R9.4.1 flowcells.
New models trained with taiyaki can be used with megalodon, but in order to obtain the highest performance the megalodon (SNP and modified base) calibration files should be reproduced for any new model.

The included model contains 5mC and 6mA capabilities.
5mC was trained only in the E. coli (CCWGG) and human (CpG) contexts while the 6mA was trained only on the E. coli (GATC) context.
Modified base detection outside of these contexts has not been tested and may produce sub-par results.
As noted above newly trained models using taiyaki can be used with megalodon, but calibration files should be reproduced for each new model.
