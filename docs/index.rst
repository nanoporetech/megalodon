*************************************
Welcome to Megalodon's documentation!
*************************************

Megalodon provides "basecalling augmentation" for raw nanopore sequencing reads, including direct, reference-guided SNP and modified base calling.

Megalodon anchors the information rich neural network basecalling output to a reference genome. Variants, either modified bases or alternative canonical bases, are then proposed and scored in order to produce highly-accurate reference anchored calls.

-------------
Prerequisites
-------------

Megalodon requires `taiyaki <https://github.com/nanoporetech/taiyaki>`_ installation for basecalling backend.
Megalodon requires only a minimal taiyaki installation via ``pip install git+https://github.com/nanoporetech/taiyaki.git``.
Full installation (via ``make install``) is not necessary for megalodon functionality.

Megalodon requires `pytorch <https://pytorch.org/>`_ to support the ``taiyaki`` basecalling backend.
For megalodon GPU support, pytorch must be installed with GPU support (and ``--devices`` to use provided at run time).
If pytorch is not installed before megalodon, pip will install the defualt pytorch (possibly CPU only).

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

Megalodon is accessed via the command line interface, ``megalodon`` command.

::

   # megalodon help (common args)
   megalodon -h
   # megalodon help (all args)
   megalodon --help-long

   # Example command calling variants and CpG methylation
   #   Compute settings: GPU devices 0 and 1 with 8 CPU cores
   megalodon raw_fast5s/ \
       --outputs basecalls mappings snps mods \
       --reference reference.fa --variant-filename variants.vcf \
       --mod-motif Z CG 0 --devices 0 1 --processes 8 --verbose-read-progress 3

This command produces the ``megalodon_results`` output directory containing basecalls, mappings, SNP and modified base results.

The majority of megalodon's functionality is accessed via the ``megalodon`` command (exemplified above), though a small number of additional scripts are found in the ``scripts`` directory of the code repository.
These include independent modified base or SNP aggregation (much faster than re-computing per-read calls), modified base result validation, and model statistic calibration.
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
