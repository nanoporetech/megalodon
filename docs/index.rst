*************************************
Welcome to Megalodon's documentation!
*************************************

Megalodon is a research command line tool to extract high accuracy modified base and sequence variant calls from raw nanopore reads by anchoring the information rich basecalling neural network output to a reference genome/transriptome.

Raw nanopore reads are processed by a single command to produce basecalls (FASTA/Q), reference mappings (SAM/BAM/CRAM), sequence variant calls (per-read and VCF) and modified base calls (per-read and bedgraph/bedmethyl/modVCF).

-------------
Prerequisites
-------------

As of version 2.0, the primary megalodon run mode requires the guppy basecaller.
See the `community page for download/installation instructions <https://community.nanoporetech.com/downloads>`_.

All other requirements are handled by ``pip`` or ``conda`` installation.
If installing from source, ``numpy`` must be installed before running installation for cython optimizations.
Required python packages are: ``cython``, ``h5py``, ``mappy``, ``numpy``, ``ont_fast5_api``, ``pyguppyclient``, ``pysam``, ``scipy``, and ``tqdm``.

.. note::

   `Taiyaki <https://github.com/nanoporetech/taiyaki>`_ is no longer required to run megalodon, but installation is required for two specific run modes:
   1) output mapped signal files (for basecall models training)
   2) runing the taiyaki basecalling backend (for neural network designs including experimental layers)

------------
Installation
------------

Megalodon is a command line tool.
``pip`` and ``conda`` are the recommended installation interfaces for megalodon.

::

   pip install megalodon
   # or
   conda install megalodon

===========
Quick Start
===========

Megalodon is accessed via the command line interface ``megalodon`` command.
The path to the ``guppy_basecall_server`` executable is required to run megalodon.
By default, megalodon assumes this path is ``./ont-guppy/bin/guppy_basecall_server``.
Use the ``--guppy-server-path`` argument to specify a different path.

::

    # megalodon help (common args)
    megalodon -h
    # megalodon help (all args)
    megalodon --help-long

    # Example command to output basecalls, mappings, variants and CpG methylation
    #   Compute settings: GPU devices 0 and 1 with 40 CPU cores
    megalodon raw_fast5s/ \
        --outputs basecalls mappings variants mods \
        --reference reference.fa --variant-filename variants.vcf.gz \
        --mod-motif Z CG 0 --devices 0 1 --processes 40 \
        --verbose-read-progress 3

This command produces the ``megalodon_results`` output directory containing all requested output files and logs.

The majority of megalodon's functionality is accessed via the ``megalodon`` command (exemplified above), though a small number of additional scripts are found in the ``scripts`` directory of the code repository.
These scripts include modified base or variant aggregation (much faster than re-computing per-read calls), modified base result validation, and model statistic calibration.
Helper scripts to perform sequence variant phasing (details here :doc:`variant_phasing`) are also included in the ``scripts`` directory of the repository.
In the future these script will move to a dedicated command line interface (likely ``megalodon_extras``).

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
   model_training
   modbase_training
