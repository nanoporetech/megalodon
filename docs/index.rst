*************************************
Welcome to Megalodon's documentation!
*************************************

Megalodon provides "basecalling augmentation" for raw nanopore sequencing reads, including direct, reference-guided SNP and modified base calling.

Megalodon anchors the information rich neural network basecalling output to a reference genome. Variants, either modified bases or alternative bases, are then proposed and scored in order to produce highly-accurate reference anchored calls.

------------
Installation
------------

.. |pypi_badge| image:: https://badge.fury.io/py/ont-megalodon.svg
    :target: https://badge.fury.io/py/ont-megalodon

Basic megalodon installation (python 3 support only)

::

    git clone https://git.oxfordnanolabs.local/algorithm/megalodon
    cd megalodon
    pip install .

See :doc:`tutorials` for common workflows.

===========
Quick Start
===========

Megalodon is accessed via the command line interface, ``megalodon`` command.

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

The majority of megalodon's functionality is accessed via the ``megalodon`` command (exemplified above), though a small number of additional scripts are found in the ``scripts`` directory of the code repository.
These including independent modified base or SNP aggregation (much faster than per-read calls), result validation, and model statistic calibration.

--------
Contents
--------

.. toctree::
   :maxdepth: 2

   algorithm_details
   common_arguments
   advanced_arguments
   computing_considerations
