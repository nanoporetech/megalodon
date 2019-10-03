.. image:: /ONT_logo.png
  :width: 800

******************

Megalodon
"""""""""

Megalodon provides "basecalling augmentation" for raw nanopore sequencing reads, including direct, reference-guided sequence variant and modified base calling.

Megalodon anchors the information rich neural network basecalling output to a reference genome.
Variants, modified bases or alternative canonical bases, are then proposed and scored in order to produce highly-accurate reference anchored modified base or sequence variant calls.

Detailed documentation for all ``megalodon`` arguments and algorithms can be found on the `megalodon documentation page <https://nanoporetech.github.io/megalodon/>`_.

Prerequisites
-------------

Megalodon requires ``numpy`` to be installed before running megalodon installation command (install with ``pip install numpy``).

Megalodon requires `taiyaki <https://github.com/nanoporetech/taiyaki>`_ installation for basecalling backend at run time.
Megalodon requires only a minimal taiyaki installation via ``pip install git+https://github.com/nanoporetech/taiyaki.git``.

Megalodon requires `pytorch <https://pytorch.org/>`_ to support the ``taiyaki`` basecalling backend.
For megalodon GPU support, pytorch must be installed with GPU support (and ``--devices`` to use provided at run time).
If pytorch is not installed before megalodon, pip will install the defualt pytorch (possibly CPU only).

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

    # Example command calling variants and CpG methylation
    #   Compute settings: GPU devices 0 and 1 with 8 CPU cores
    megalodon raw_fast5s/ \
        --outputs basecalls mappings variants mods \
        --reference reference.fa --variant-filename variants.vcf.gz \
        --mod-motif Z CG 0 --devices 0 1 --processes 8 --verbose-read-progress 3

This command produces the ``megalodon_results`` output directory containing basecalls, mappings, sequence variants and modified base results.
The format for each output is described below.

.. note::

   The default basecalling model (used in this example command) is for R9.4.1, MinION/GridION reads.
   This is equivalent to the guppy/MinKNOW high accuracy model in terms of basecall accuracy and run speed (though using taiyaki backend is slower than guppy).
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

    - If not indexed, indexing will be performed
  - Megalodon currently requires a set of candidate variants for ``--outputs variants``.
  - Only small indels (default less than ``50`` bases) are tested by default.

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
  - A tab-separated mapping text summary is produced including per-read alignment statistics

    - ``percent_identity`` is defined as ``num_matched_bases`` / ``num_align_bases``
- Modified Base Calls

  - Per-read modified base calls

    - Per-read SQL DB containing scores at each tested reference location

      - Contains an indexed table with per-read, per-position, modified base scores, as well as auxiliary tables with read, modification type and reference position information.
    - Tab-delimited output can be produced by adding the ``--write-mods-text`` flag
  - Aggregated calls

    - Aggregated calls are output in either bedMethyl format (default; one file per modified base), a VCF variant format (including all modified bases) or wiggle format (one file per modified base/strand combination).
- Sequence Variant Calls

  - Per-read Variant Calls

    - SQL DB containing scores for each tested variant

      - Contains a single ``variants`` table indexed by reference position
    - Tab-delimited output can be produced by adding the ``--write-variants-text`` flag
  - Aggregated calls

    - Format: VCF
    - VCF sample field contains ``GT``, ``GQ``, ``DP``, ``GL``, and ``PL`` attributes
    - Default run mode is diploid. To run in haploid mode, set ``--haploid`` flag.

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
These include enumerating read ids and files, collecting and reporting progress information and getting data from read processing queues and writing outputs (basecalls, mappings, sequence variants and modified bases).

Model Compatibility
-------------------

The model and calibration files included with megalodon are applicable only to MinION or GridION R9.4.1 flowcells.
New models trained with taiyaki can be used with megalodon, but in order to obtain the highest performance the megalodon (variant and modified base) calibration files should be reproduced for any new model (TODO provide walkthrough).

The included model contains 5mC and 6mA capabilities.
5mC was trained only in the human (CpG) and E. coli (CCWGG) contexts while the 6mA was trained only on the E. coli (GATC) context.
Modified base detection outside of these contexts has not been tested and may produce sub-optimal results.
As noted above newly trained models using taiyaki can be used with megalodon, but calibration files should be reproduced for each new model.

RNA
---

Megalodon does not currently support direct RNA processing.
This feature is currently under development.

Licence and Copyright
---------------------

|copy| 2019 Oxford Nanopore Technologies Ltd.

.. |copy| unicode:: 0xA9 .. copyright sign

Megalodon is distributed under the terms of the Oxford Nanopore
Technologies, Ltd.  Public License, v. 1.0.  If a copy of the License
was not distributed with this file, You can obtain one at
http://nanoporetech.com
