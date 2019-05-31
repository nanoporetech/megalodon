.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

Megalodon provides "basecalling augmentation" for raw nanopore sequencing reads, including direct, reference-guided SNP and modified base calling.

Megalodon anchors the information rich neural network basecalling output to a reference genome.
Variants, modified bases or alternative canonical bases, are then proposed and scored in order to produce highly-accurate reference anchored modified base or SNP calls.

Detailed documentation for all ``megalodon`` arguments and algorithms can be found on the `megalodon documentation page <https://nanoporetech.github.io/megalodon/>`_.

Prerequisites
-------------

Megalodon requires ``numpy`` and ``cython`` to be installed before running megalodon installation command.

Megalodon requires `taiyaki <https://github.com/nanoporetech/taiyaki>`_ installation for basecalling backend at run time.
Megalodon requires only a minimal taiyaki installation via ``pip install git+git://github.com/nanoporetech/taiyaki.git``.
Full installation via ``git clone https://github.com/nanoporetech/taiyaki && cd taiyaki && make install`` is not necessary for megalodon functionality.

Megalodon requires `pytorch <https://pytorch.org/>`_ to support the ``taiyaki`` basecalling backend.
For megalodon GPU support, pytorch must be installed with GPU support (and ``--devices`` to use provided at run time).
If pytorch is not installed before megalodon, pip will install the defualt pytorch (possibly CPU only).

Installation
------------

::

   pip install numpy cython
   git clone https://github.com/nanoporetech/megalodon
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

    # Example command calling variants and CpG methylation
    #   Compute settings: GPU devices 0 and 1 with 8 CPU cores
    megalodon raw_fast5s/ \
        --outputs basecalls mappings snps mods \
        --reference reference.fa --variant-filename variants.vcf.gz \
        --mod-motif Z CG 0 --devices 0 1 --processes 8 --verbose-read-progress 3

This command produces the ``megalodon_results`` output directory containing basecalls, mappings, SNP and modified base results.
The format for each output is described below.

.. note::

   The default basecalling model (used in this example command) is for R9.4.1, MinION/GridION reads.
   This is equivalent to the guppy/MinKNOW high accuracy model in terms of basecall accuracy and run speed (though using taiyaki backend is slower than guppy).
   This model contains 5mC (encoded as a ``Z`` base) and 6mA (encoded as a ``Y`` base) trained in biological contexts only (5mC in human CpG and E. coli CCWGG and 6mA in E. coli GATC).

Inputs
------

- Raw reads

  - Directory containing raw read FAST5 files
  - By default the directory will be searched recursively for read files (ending in ``.fast5``)
- Reference

  - Genome or transcriptome sequence reference file in FASTA format
- Variants File

  - Format: indexed VCF or BCF
  - Optional, but required for SNP calling
  - Megalodon currently requires a set of candidate variants in order to call SNPs.
  - Only small indels (default less than ``50`` bases) are included in testing.

    - Specify the ``--max-indel-size`` argument to process larger indels
    - The ``--variant-context-bases`` argument may need to be increased for larger indels.

Outputs
-------

- Basecalls

  - Format: FASTA

    - FASTQ format output is not currently available
  - Basecall-anchored modified base scores are also available (via HDF5 output)
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

    - Phased VCF output
    - Phased read calls
    - Improved phase-aware per-read sequence variants calls

Computing
---------

Megalodon processes reads from a queue using a pool of workers.
The number of workers is set using the ``--processes`` argument.
Each process is linked to a taiyaki basecalling backend and a separate thread for reference mapping.
The threaded mapping interface allows megalodon to load the reference (via ``mappy``) into shared memory.

In order to use GPU resources the ``--devices`` argument can be set.
If ``--devices`` is set, the taiyaki backends will be distribured evenly over the specified ``--devices``.
In order to control the GPU memory usage, the ``--max-concurrent-chunks`` argument allows a user to restrict the maximum number of chunks to process concurrently (per ``--process``).
Note that the model parameters must (currently) be loaded into each GPU process and thus limits the number of GPU processes that can be spawned per GPU.

The ``--chunk-size`` and ``--chunk-overlap`` arguments allow users to specify read chunking, but signal normalization is always carried out over the entire read.

A number of helper processes will be spawned in order to perform more minor tasks, which should take minimal compute resources.
These include enumerating read ids and files, collecting and reporting progress information and getting data from read processing queues and writing outputs (basecalls, mappings, SNPs and modified bases).

Compatibility
-------------

The model and calibration files included with megalodon are applicable only to MinION or GridION R9.4.1 flowcells.
New models trained with taiyaki can be used with megalodon, but in order to obtain the highest performance the megalodon (SNP and modified base) calibration files should be reproduced for any new model (TODO provide walkthrough).

The included model contains 5mC and 6mA capabilities.
5mC was trained only in the human (CpG) and E. coli (CCWGG) contexts while the 6mA was trained only on the E. coli (GATC) context.
Modified base detection outside of these contexts has not been tested and may produce sub-optimal results.
As noted above newly trained models using taiyaki can be used with megalodon, but calibration files should be reproduced for each new model.


Licence and Copyright
---------------------

|copy| 2019 Oxford Nanopore Technologies Ltd.

.. |copy| unicode:: 0xA9 .. copyright sign

Megalodon is distributed under the terms of the Oxford Nanopore
Technologies, Ltd.  Public License, v. 1.0.  If a copy of the License
was not distributed with this file, You can obtain one at
http://nanoporetech.com
