*************************************
Welcome to Megalodon's documentation!
*************************************

Megalodon is a research tool for per-read and aggregated modified base and sequence variant calling by anchoring the information rich basecalling neural network output to a reference genome/transriptome.
Megalodon takes raw nanopore reads as input and produces multiple outputs, primarily including basecalls (FASTA), reference mappings (SAM/BAM/CRAM), sequence variant calls (per-read and VCF) and modified base calls (per-read and bedgraph/bedmethyl/VCF).

-------------
Prerequisites
-------------

Megalodon requires ``numpy`` to be installed before running megalodon installation command (install with ``pip install numpy``).

As of version 1.0.0, megalodon can be run directly on FAST5 files output by guppy (requires ``--fast5_out`` and ``--post_out`` flags).
Using guppy FAST5 input, megalodon requires no other prerequisites (aside from those automatically installed by ``pip``).

In order to run megalodon directly from raw (un-basecalled) FAST5 files, `taiyaki <https://github.com/nanoporetech/taiyaki>`_ is currently required.
Megalodon requires only a minimal taiyaki installation via ``pip install git+https://github.com/nanoporetech/taiyaki.git``.
Full ``taiyaki`` installation (via ``make install``) is not necessary for megalodon functionality, but makes GPU configuration (via ``pytorch``) more painless.

.. note::

   The two currently implemented basecalling backends (``guppy`` and ``taiyaki``) have a number of advantages and disadvantages.
   ``Guppy`` is computationally faster but requires disk storage of the neural network outputs (~10 times increase over raw reads).
   The ``guppy`` backend should be chosen if speed is important and a sufficient amount of temporary fast disk space is available.
   The ``taiyaki`` basecalling backend is computationally slower than ``guppy`` (usually 5-10 times) but requires storing only the minimal outputs on disk.
   For very advanced users, the ``taiyaki`` backend is also extensible to any ``taiyaki`` flip-flop model architecture while ``guppy`` only supports a set of optimized basecalling model architectures.

------------
Installation
------------

Basic megalodon installation (python 3 support only)

::

   pip install numpy
   pip install git+https://git.oxfordnanolabs.local/algorithm/megalodon.git

===========
Quick Start
===========

Megalodon is accessed via the command line interface ``megalodon`` command.

::

    # megalodon help (common args)
    megalodon -h
    # megalodon help (all args)
    megalodon --help-long

    # Example commands calling variants and CpG methylation
    #   Compute settings: GPU devices 0 and 1 with 10 CPU cores

    # Using taiyaki backend
    megalodon raw_fast5s/ \
        --load-default-model --outputs basecalls mappings variants mods \
        --reference reference.fa --variant-filename variants.vcf.gz \
        --mod-motif Z CG 0 --devices 0 1 --processes 10 --verbose-read-progress 3

    # Using guppy backend
    ont-guppy/bin/guppy_basecaller \
        -i raw_fast5s/ -s basecalled_fast5s/ \
        -c ont-guppy/data/dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg \
       --post_out --fast5_out --device cuda:0 --num_callers 10
    megalodon basecalled_fast5s/ \
        --outputs basecalls mappings variants mods \
        --reference reference.fa --variant-filename variants.vcf.gz \
        --mod-motif Z CG 0 --devices 0 1 --processes 10 --verbose-read-progress 3

This command produces the ``megalodon_results`` output directory containing basecalls, mappings, sequence variant and modified base results.

The majority of megalodon's functionality is accessed via the ``megalodon`` command (exemplified above), though a small number of additional scripts are found in the ``scripts`` directory of the code repository.
These include modified base or variant aggregation (much faster than re-computing per-read calls), modified base result validation, and model statistic calibration.
Helper scripts to perform sequence variant phasing (details here :doc:`variant_phasing`) are also included in the ``scripts`` directory of the repository.

--------
Contents
--------

.. toctree::
   :maxdepth: 2

   algorithm_details
   common_arguments
   advanced_arguments
   computing_considerations
   variant_phasing
   file_formats
