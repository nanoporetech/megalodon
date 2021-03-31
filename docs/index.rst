*************************************
Welcome to Megalodon's documentation!
*************************************

Megalodon is a research command line tool to extract high accuracy modified base and sequence variant calls from raw nanopore reads by anchoring the information rich basecalling neural network output to a reference genome/transcriptome.

Raw nanopore reads are processed by a single command to produce basecalls (FASTA/Q), reference mappings (SAM/BAM/CRAM), modified base calls (per-read and aggregated per-reference site), sequence variant calls (per-read and aggregated per-reference site) and more.

-------------
Prerequisites
-------------

The primary Megalodon run mode requires the Guppy basecaller (version >= 4.0).
See the `community page for download/installation instructions [login required] <https://community.nanoporetech.com/downloads>`_.

Megalodon is a python-based command line software package.
Given a python (version >= 3.5) installation, all other requirements are handled by ``pip`` or ``conda``.

..

   `Taiyaki <https://github.com/nanoporetech/taiyaki>`_ is no longer required to run Megalodon, but installation is required for two specific run modes:

   1) output mapped signal files (for basecall model training)

   2) running the Taiyaki basecalling backend (for neural network designs including experimental layers)

------------
Installation
------------

``pip`` is recommended for Megalodon installation.

::

   pip install megalodon

``conda`` installation is available, but not fully supported.
``ont_pyguppy_client_lib`` is not available on conda and thus must be installed with ``pip``.

::

   conda install megalodon
   pip install ont_pyguppy_client_lib

To install from github source for development, the following commands can be run.

::

   git clone https://github.com/nanoporetech/megalodon
   pip install -e megalodon/

It is recommended that Megalodon be installed in a control compute environment.
See `the python documentation for preparing virtual environments <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_

===========
Quick Start
===========

Megalodon must obtain the intermediate output from the basecall neural network.
Guppy (production nanopore basecalling software) is the recommended backend to obtain this output from raw nanopore signal (from FAST5 files).
Nanopore basecalling is compute intensive and thus it is highly recommended that GPU resources are specified (``--devices``) for optimal Megalodon performance.

Megalodon is accessed via the command line interface ``megalodon`` command.

::

    # megalodon help (common args)
    megalodon -h
    # megalodon help (advanced args)
    megalodon --help-long

    # Example command to output basecalls, mappings, and 5mC CpG methylation in both per-read (``mod_mappings``) and aggregated (``mods``) formats
    #   Compute settings: GPU devices 0 and 1 with 40 CPU cores
    megalodon \
        raw_fast5s/ \
        --outputs basecalls mappings mod_mappings mods \
        --reference reference.fa --mod-motif m CG 0 \
        --devices 0 1 --processes 40

This command produces the ``megalodon_results`` output directory containing all requested output files and logs.
The format for common outputs is described briefly below and in more detail in the `full documentation <https://nanoporetech.github.io/megalodon/>`_

The above command uses the modified base model included in Guppy (more details below `Guppy Models and Parameters`_).
As of the ``2.3.0`` megalodon release (March 2020) the models included with Guppy provide the most accurate modified basecalling models.
As more accurate basecalling models are trained, they are first released into the `Rerio repository for research models <https://github.com/nanoporetech/rerio>`_.
Once training pipelines are more thoroughly standardized and tested models will be transferred into Guppy.
The code below shows how to obtain and run the R9.4.1, MinION/GridION, 5mC CpG model from Rerio.
Note that this is the same model now included in Guppy 4.5.0+.

::

    # Obtain and run R9.4.1, MinION, 5mC CpG model from Rerio
    git clone https://github.com/nanoporetech/rerio
    rerio/download_model.py rerio/basecall_models/res_dna_r941_min_modbases_5mC_CpG_v001
    megalodon \
        raw_fast5s/ \
        --guppy-params "-d ./rerio/basecall_models/" \
        --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg \
        --outputs basecalls mappings mod_mappings mods \
        --reference reference.fa --mod-motif m CG 0 \
        --devices 0 1 --processes 40

..

    The path to the ``guppy_basecall_server`` executable is required to run Megalodon.
    By default, Megalodon assumes Guppy (Linux GPU) is installed in the current working directory (i.e. ``./ont-guppy/bin/guppy_basecall_server``).
    Use the ``--guppy-server-path`` argument to specify a different path.

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
   extras_aggregate
   extras_calibrate
   extras_merge
   extras_modified_bases
   extras_phase_variants
   extras_per_read_text
   extras_validate
   extras_variants
