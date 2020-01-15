.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

Megalodon is a research tool for per-read and aggregated modified base and sequence variant calling by anchoring the information rich basecalling neural network output to a reference genome/transriptome.
Megalodon takes raw nanopore reads as input and produces multiple outputs, primarily including basecalls (FASTA), reference mappings (SAM/BAM/CRAM), sequence variant calls (per-read and VCF) and modified base calls (per-read and bedgraph/bedmethyl/VCF).

Detailed documentation for all ``megalodon`` arguments and algorithms can be found on the `megalodon documentation page <https://nanoporetech.github.io/megalodon/>`_.

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

Installation
------------

::

   pip install numpy cython
   pip install git+https://github.com/nanoporetech/megalodon.git

Getting Started
---------------

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

This command produces the ``megalodon_results`` output directory containing basecalls, mappings, sequence variants and modified base results.
The format for each output is described below and in more detail on the `documentation page <https://nanoporetech.github.io/megalodon/>`_

.. note::

   The default basecalling model (used in this example command) is for R9.4.1, MinION/GridION reads.
   This model contains modified bases 5mC (encoded as a ``Z`` base) and 6mA (encoded as a ``Y`` base) trained in biological contexts only (5mC in human CpG and E. coli CCWGG and 6mA in E. coli GATC).

Inputs
------

- Raw reads

  - Directory containing raw read FAST5 files
  - By default the directory will be searched recursively for read files (ending in ``.fast5``)
- Reference

  - Genome or transcriptome sequence reference file in FASTA format (minimap2 index may be provided instead)
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

  - Format: FASTA

    - FASTQ format output is not currently available
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

High Quality Phased Variant Calls
---------------------------------

In order to obtain the highest quality diploid sequence variant calls the full variant phasing pipeline employing whatshap should be applied.
This pipeline is described in detail on the `full documentation page <https://nanoporetech.github.io/megalodon/variant_phasing.html>`_.
The default diploid variant settings are optimized for this full phasing pipeline and not for direct validation at this point.
This includes overcalling heterozygous sites in order to accurately phase as many variants as possible via whatshap.
Thus validation of the direct diploid variant calling results will likely show overcalling of heterozygous sites.

Taiyaki-backend Computing
-------------------------

Megalodon processes reads from a queue using a pool of workers.
The number of workers is set using the ``--processes`` argument.
Each process is linked to a taiyaki basecalling backend and a separate thread for reference mapping.
The threaded mapping interface allows megalodon to load the reference into shared memory (via `mappy <https://github.com/lh3/minimap2/tree/master/python>`_).

In order to use GPU resources set the ``--devices`` argument.
If ``--devices`` is set, the taiyaki backends will be distribured evenly over the specified ``--devices``.
In order to control the GPU memory usage, the ``--max-concurrent-chunks`` argument allows a user to restrict the maximum number of chunks to process concurrently (per ``--process``).
Note that the model parameters must be loaded into each GPU process (default high accuracy model consumes ~1GB of GPU memory) and thus limits the number of processes that can be spawned per GPU device.

The ``--chunk-size`` and ``--chunk-overlap`` arguments allow users to specify read chunking, but signal normalization is always carried out over the entire read.

A number of helper processes will be spawned in order to perform more minor tasks, which should take minimal compute resources.
These include enumerating read ids and files, collecting and reporting progress information and getting data from read processing queues and writing outputs (basecalls, mappings, sequence variants and modified bases).

Disk Performance Considerations
*******************************

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

The model and calibration files included with megalodon are applicable only to MinION or GridION R9.4.1 flowcells.
New models trained with taiyaki can be used with megalodon, but in order to obtain the highest performance the megalodon (variant and modified base) calibration files should be reproduced for any new model.

The default model included with megalodon provides 5mC and 6mA detection.
5mC was trained only in the human (CpG) and E. coli (CCWGG) contexts while the 6mA was trained only on the E. coli (GATC) context.
Modified base detection outside of these contexts has not been tested and may produce sub-optimal results.
As noted above newly trained models using taiyaki can be used with megalodon, but calibration files should be reproduced for each new model.

RNA
---

Megalodon provides experimental support for direct RNA processing.
This support can be accessed within the ``rna`` code branch (access via ``git clone --branch rna https://github.com/nanoporetech/megalodon``).

Licence and Copyright
---------------------

|copy| 2019 Oxford Nanopore Technologies Ltd.

.. |copy| unicode:: 0xA9 .. copyright sign

Megalodon is distributed under the terms of the Oxford Nanopore
Technologies, Ltd.  Public License, v. 1.0.  If a copy of the License
was not distributed with this file, You can obtain one at
http://nanoporetech.com
