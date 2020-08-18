.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

.. image:: https://img.shields.io/pypi/v/megalodon   :alt: PyPI
.. image:: https://img.shields.io/pypi/dm/megalodon   :alt: PyPI - Downloads

.. image:: https://img.shields.io/conda/vn/bioconda/megalodon   :alt: Conda
.. image:: https://img.shields.io/conda/dn/bioconda/megalodon   :alt: Conda - Downloads

Megalodon is a research command line tool to extract high accuracy modified base and sequence variant calls from raw nanopore reads by anchoring the information rich basecalling neural network output to a reference genome/transcriptome.

Raw nanopore reads are processed by a single command to produce basecalls (FASTA/Q), reference mappings (SAM/BAM/CRAM), modified base calls (per-read and bedgraph/bedmethyl/modVCF), sequence variant calls (per-read and VCF) and more.

Detailed documentation for all ``megalodon`` commands and algorithms can be found on the `Megalodon documentation page <https://nanoporetech.github.io/megalodon/>`_.

Prerequisites
-------------

As of version 2.0, the primary Megalodon run mode requires the Guppy basecaller and as of version 2.2 Guppy version >= 4.0 is required.
See the `community page for download/installation instructions [login required] <https://community.nanoporetech.com/downloads>`_.

Megalodon is a python-based command line software package.
Given a python (version >= 3.5) installation, all other requirements are handled by ``pip`` or ``conda``.

..

   `Taiyaki <https://github.com/nanoporetech/taiyaki>`_ is no longer required to run megalodon, but installation is required for two specific run modes:
   1) output mapped signal files (for basecall model training)
   2) running the Taiyaki basecalling backend (for neural network designs including experimental layers)

Installation
------------

``pip`` and ``conda`` are the recommended installation interfaces for Megalodon.
``ont_pyguppy_client_lib`` is not available on conda and thus must be installed with ``pip``.

::

   pip install megalodon
   # or
   conda install megalodon
   pip install ont_pyguppy_client_lib

To install from github source for development, the following commands can be run.
``numpy`` must be installed before running installation for cython optimizations.

::

   git clone https://github.com/nanoporetech/megalodon
   pip install numpy cython
   pip install -e megalodon/

Getting Started
---------------

Megalodon is accessed via the command line interface ``megalodon`` command.

::

    # megalodon help (common args)
    megalodon -h
    # megalodon help (all args)
    megalodon --help-long

    # Example command to output basecalls, mappings, and CpG methylation in both per-read (``mod_mappings``) and aggregated (``mods``) formats
    #   Compute settings: GPU devices 0 and 1 with 40 CPU cores
    # For highest accuracy methylation calls see research models in Rerio: https://github.com/nanoporetech/rerio
    megalodon raw_fast5s/ \
        --outputs basecalls mappings mod_mappings mods \
        --reference reference.fa --mod-motif Z CG 0 \
        --devices 0 1 --processes 40

This command produces the ``megalodon_results`` output directory containing all requested output files and logs.
The format for common outputs is described briefly below and in more detail in the `full documentation <https://nanoporetech.github.io/megalodon/>`_

..

    The path to the ``guppy_basecall_server`` executable is required to run Megalodon.
    By default, Megalodon assumes Guppy (Linux GPU) is installed in the current working directory (i.e. ``./ont-guppy/bin/guppy_basecall_server``).
    Use the ``--guppy-server-path`` argument to specify a different path.

Inputs
------

- Raw reads

  - Directory containing raw read FAST5 files (sub-directories recursively searched)
- Reference

  - Genome or transcriptome sequence reference (FASTA or minimap2 index)
- Variants File

  - Megalodon requires a set of candidate variants for ``--outputs variants`` (provide via ``--variant-filename`` argument; VCF or BCF).

Outputs
-------

All Megalodon outputs are written into the directory specified with the ``--output-directory`` option with standard file names and extensions.

- Basecalls

  - Format: FASTQ (default) or FASTA
  - Basecall-anchored modified base scores are also available in hts-spec BAM format tags (``--outputs mod_basecalls``).
- Mappings

  - Format: SAM, BAM (default), or CRAM
  - A tab-separated mapping text summary is also produced including per-read alignment statistics.
- Modified Base Calls

  - The basecalling model specifies the modified bases capable of being output. See ``megalodon_extras modified_bases describe_alphabet``.
  - Per-read modified base calls

    - SQL DB containing per-read modified base scores at each covered reference location
    - Reference-anchored per-read modified base calls is BAM format via the ``Mm`` and ``Ml`` tags (see `hts-spec specifications here <https://github.com/samtools/hts-specs/pull/418>`_).
  - Aggregated calls

    - Format: bedgraph, bedmethyl (default), and/or modVCF
  - In order to restrict modified base calls to a specific motif(s) specify the ``--mod-motif`` argument. For example, to restrict calls to CpG sites specify ``--mod-motif Z CG 0``.
- Sequence Variant Calls

  - Per-read Variant Calls

    - SQL DB containing per-read variant scores for each covered variant
  - Aggregated calls

    - Format: VCF
    - Default run mode is diploid. To run in haploid mode, set ``--haploid`` flag.
    - For best results on a diploid genome see the variant phasing workflow on the `full documentation page <https://nanoporetech.github.io/megalodon/variant_phasing.html>`_.

Live Processing
---------------

As of version 2.2, Megalodon now supports live run processing.
Activate live processing mode by simply adding the ``--live-processing`` argument and specifying the MinKNOW output directory as the Megalodon FAST5 input directory.
Megalodon will continue to search for FAST5s until the ``final_summary*`` file is created by MinKNOW, indicating data production has completed.

Guppy Models and Parameters
---------------------------

As of version 2.2, Megalodon requires Guppy version >= 4.0.

The Guppy model defines the modified bases capable of being output by Megalodon.
Basecalling models must be trained to specifically detect a type or types of modified bases.
See the `Megalodon documentation here <https://nanoporetech.github.io/megalodon/modbase_training.html>`_ for instructions to construct modified base training data and train a new modified base model.

By default, Megalodon uses the ``dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg`` Guppy config.
This config is compatible with DNA, R9.4.1, MinION/GridION reads and allows output of 5mC and 6mA calls in biological contexts (CpG, dcm and dam sites).
Use the ``--guppy-config`` option to specify a different guppy model config.
The appropriate `Rerio model <https://github.com/nanoporetech/rerio>`_ is recommended for the highest accuracy modified base calls.

All configs can be used to output ``basecalls`` and ``mappings`` (as well as ``signal_mappings`` and ``per_read_refs`` for `basecall training <https://nanoporetech.github.io/megalodon/model_training.html>`).
Modified base and sequence variant outputs require Megalodon calibration files.
To list configs with default calibration files, run ``megalodon --list-supported-guppy-configs``.
See `calibration documentation here <https://nanoporetech.github.io/megalodon/extras_calibrate.html>`_ for details on Megalodon model calibration.

Only flip-flop configs/models are currently supported by Megalodon (this excludes k-mer based and RLE model types).

In addition to the ``--guppy-config`` and ``--guppy-server-path`` options, a number of additional arguments control the behavior of the guppy backend.
The ``--guppy-params`` argument will pass arguments directly to the ``guppy_basecall_server`` initialization call.
For example to optimize GPU usage, the following option might be specified: ``--guppy-params "--num_callers 5 --ipc_threads 6"``

Finally the ``--guppy-timeout`` arguments ensures that a run will not stall on a small number of reads taking a very long time (default 5 seconds).

High Quality Phased Variant Calls
---------------------------------

In order to obtain the highest quality diploid sequence variant calls, the full variant phasing pipeline employing ``whatshap`` should be applied.
This pipeline is described in detail on the `full documentation page <https://nanoporetech.github.io/megalodon/variant_phasing.html>`_.
The default diploid variant settings are optimized for the full phasing pipeline and not the highest quality diploid calls directly from a single Megalodon call.

High-Density Variants
---------------------

When running Megalodon with a high density of variants (more than 1 variant per 100 reference bases), certain steps can be taken to increase performance.
See `variant atomize documentation <https://nanoporetech.github.io/megalodon/extras_variants.html#megalodon-extras-variants-atomize>`_ for further details.

Disk Performance Considerations
-------------------------------

Per-read modified base and variant statistics are stored in an on-disk sqlite database.
As of version 2.0, the status of output queues and as of version 2.2 the extract signal input queue are displayed by default.
If the ``extract_signal`` input queue is often empty, Megalodon is waiting on reading raw signal from FAST5 files.
If any output status bars indicate a full queue, Megalodon will stall waiting on that process to write data to disk.
Moving the input data directory or  ``--output-directory`` accordingly to a location with faster disk I/O performance should improve performance.

RNA
---

Megalodon supports processing direct RNA nanopore data.
In order to process an RNA sample specify the ``--rna`` flag as well as an RNA model using the ``--guppy-config`` argument.

Megalodon performs mapping using the standard minimap2 option, ``map-ont``, and not the ``splice`` option, so a transcriptome reference must be provided.
The Megalodon code supports RNA modified base detection, but currently no RNA modified base basecalling models are released.

.. note::

   Megalodon does not currently perform checking that a set of reads agree with the provided model or options specified (e.g. ``--rna``).
   Users should take care to ensure that the correct options are specified for each sample processed.

License and Copyright
---------------------

|copy| 2019-20 Oxford Nanopore Technologies Ltd.

.. |copy| unicode:: 0xA9 .. copyright sign

Megalodon is distributed under the terms of the Oxford Nanopore
Technologies, Ltd.  Public License, v. 1.0.  If a copy of the License
was not distributed with this file, You can obtain one at
http://nanoporetech.com
