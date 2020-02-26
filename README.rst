.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

Megalodon is a research tool to extract high accuracy modified base and sequence variant calls from raw nanopore reads by anchoring the information rich basecalling neural network output to a reference genome/transriptome.
Raw nanopore reads are processed by a single command to produce basecalls (FASTA/Q), reference mappings (SAM/BAM/CRAM), sequence variant calls (per-read and VCF) and modified base calls (per-read and bedgraph/bedmethyl/modVCF).

Detailed documentation for all ``megalodon`` arguments and algorithms can be found on the `megalodon documentation page <https://nanoporetech.github.io/megalodon/>`_.

Prerequisites
-------------

As of version 2.0, the primary megalodon run mode requires the guppy basecaller.
See the `community page for download/installation instructions <https://community.nanoporetech.com/downloads>`_.

All other requirements are handled by ``pip`` or ``conda`` installation.
If installing from source, ``numpy`` must be installed before running installation for cython optimizations.

.. note::

   `Taiyaki <https://github.com/nanoporetech/taiyaki>`_ installation is only required to output mapped signal files (for basecall models training) or to run the taiyaki basecalling backend (for neural network designs including experimental layers).

Installation
------------

Megalodon is a command line tool.
``pip`` and ``conda`` are the recommended installation interfaces for megalodon.

::

   pip install megalodon
   # or
   conda install megalodon

Getting Started
---------------

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
The format for each output is described briefly below and in more detail in the `full documentation <https://nanoporetech.github.io/megalodon/>`_

Inputs
------

- Raw reads

  - Directory containing raw read FAST5 files
  - By default, the directory will be searched recursively for read files (ending in ``.fast5``)
- Reference

  - Genome or transcriptome sequence reference (FASTA or minimap2 index)
- Variants File

  - Format: VCF or BCF
  - Megalodon currently requires a set of candidate variants for ``--outputs variants`` (provide via ``--variant-filename`` argument).

    - If not indexed, indexing will be performed
  - Only small indels (default less than ``50`` bases) are tested by default.

    - Specify the ``--max-indel-size`` argument to process larger indels
    - The ``--variant-context-bases`` argument may need to be increased for larger indels.

Outputs
-------

All megalodon outputs are output into the directory specified with the ``--output-directory`` option with standard file names and extensions.

- Basecalls

  - Format: FASTA/Q
  - Basecall-anchored modified base scores are also available (via HDF5 output)
- Mappings

  - Format: SAM, BAM (default), or CRAM
  - A tab-separated mapping text summary is produced including per-read alignment statistics

    - ``percent_identity`` is defined as ``num_matched_bases`` / ``num_align_bases``
- Modified Base Calls

  - In order to restrict modified base calls to a particular motif specify the ``--mod-motif`` along with the modified base, canonical motif and relative modified base position within the motif. For example in order to output only CpG methylation specify ``--mod-motif Z CG 0``.
  - Per-read modified base calls

    - Per-read SQL DB containing scores at each tested reference location

      - Contains an indexed table with per-read, per-position, modified base scores, as well as auxiliary tables with read, modification type and reference position information.
    - Tab-delimited output can be produced by adding the ``--write-mods-text`` flag or produced after a run using the ``megalodon/scripts/write_per_read_modified_base_text.py`` script.
  - Aggregated calls

    - Aggregated calls are output in either bedMethyl format (default; one file per modified base), a VCF variant format (including all modified bases) or wiggle format (one file per modified base/strand combination).
- Sequence Variant Calls

  - Per-read Variant Calls

    - SQL DB containing scores for each tested variant

      - Contains a single ``variants`` table indexed by reference position
    - Tab-delimited output can be produced by adding the ``--write-variants-text`` flag or produced after a run using the ``megalodon/scripts/write_per_read_sequence_variant_text.py`` script.
  - Aggregated calls

    - Format: VCF
    - VCF sample field contains ``GT``, ``GQ``, ``DP``, ``GL``, and ``PL`` attributes
    - Default run mode is diploid. To run in haploid mode, set ``--haploid`` flag.
    - For best results on a diploid genome see the variant phasing workflow on the `full documentation page <https://nanoporetech.github.io/megalodon/variant_phasing.html>`_.

Guppy Models and Parameters
---------------------------

By default, megalodon uses the ``dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg`` guppy config.
This config is compatible with DNA, R9.4.1, MinION/GridION reads and allows output of 5mC and 6mA calls in biological contexts (CpG, dcm and dam sites).
Use the ``--guppy-config`` option to specify a different guppy model config.

All configs can be used to output basecalls and mappings (as well as signal mapping files for basecall training; see ``--output-signal-mappings`` in ``--help-long``).
Modified base and sequence variant outputs require megalodon calibration files.
To list configs with default calibration files, run ``megalodon --list-supported-guppy-configs``.

Only flip-flop configs/models are currently supported by megalodon (this excludes k-mer based and RLE model types).

In addition to the ``--guppy-config`` and ``--guppy-server-path`` options, a number of additional arguments control the behavior of the guppy backend.
An alternative server port can be specified with the ``--guppy-server-port`` argument (useful when multiple megalodon/guppy_server are active on the same machine).
The ``--guppy-params`` argument will pass arguments directly to the ``guppy_basecall_server`` initialization call.
These arguments must be valid arguments for the provided guppy server executable.
For example to optimize GPU usage for an nvidia V100 GPU, the following option might be specified: ``--guppy-parmas "--num_callers 5 --ipc_threads 6"``

Finally the ``--guppy-timeout`` arguments ensures that a run will not stall on a small number of reads taking a very long time (default 5 seconds).

High Quality Phased Variant Calls
---------------------------------

In order to obtain the highest quality diploid sequence variant calls the full variant phasing pipeline employing whatshap should be applied.
This pipeline is described in detail on the `full documentation page <https://nanoporetech.github.io/megalodon/variant_phasing.html>`_.
The default diploid variant settings are optimized for this full phasing pipeline and not for direct validation at this point.
This includes overcalling heterozygous sites in order to accurately phase as many variants as possible via whatshap.
Thus validation of the direct diploid variant calling results will likely show overcalling of heterozygous sites.

Disk Performance Considerations
*******************************

TODO add note about queue status bars.

Within megalodon, per-read modified base and variant statistics are stored in an on-disk sqlite database.
During read processing per-read, per-site statistics are funneled through a single thread to handle the database input.
If the requested compute resources are not being utililized to their fullest extent during read processing slow disk write is the most likely bottleneck.
Moving the database, stored within the directory specified with the ``--output-directory`` argument, to a location with faster disk I/O performance should imporove performance.

For the aggregation stage of processing the disk read speed has a magnified effect.
During aggregation binary searches for results grouped per-site must be performed over the on-disk database.
While database optimization to reduce the disk reads has been implemented the performance for data extraction can be extremely slow for large runs.
Moving the database location from a remote or network file system to a local fast (SSD) disk can increase compute efficiency as much as 100X-1000X.

Model Compatibility
-------------------

TODO update for guppy calibration files

RNA
---

Megalodon provides experimental support for direct RNA processing.
This support can be accessed within the ``rna`` code branch (access via ``git clone --branch rna https://github.com/nanoporetech/megalodon``).

Licence and Copyright
---------------------

|copy| 2019-20 Oxford Nanopore Technologies Ltd.

.. |copy| unicode:: 0xA9 .. copyright sign

Megalodon is distributed under the terms of the Oxford Nanopore
Technologies, Ltd.  Public License, v. 1.0.  If a copy of the License
was not distributed with this file, You can obtain one at
http://nanoporetech.com
