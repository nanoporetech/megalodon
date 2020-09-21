****************
Common Arguments
****************

-----------------
Required Argument
-----------------

- ``fast5s_dir``

  - Path to directory containing raw FAST5-format nanopore reads.
  - Both single and multi FAST5 formats are supported.
  - Default searches recursively for fast5 read files. To search only one-level specify `--not-recursive`.

----------------------
Guppy Backend Argument
----------------------

- ``--guppy-config``

  - Guppy config.
  - Default: ``dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg``

- ``--guppy-server-path``

  - Path to guppy server executable.
  - Default: ``./ont-guppy/bin/guppy_basecall_server``

----------------
Output Arguments
----------------

- ``--live-processing``

  - As of version 2.2, Megalodon now supports live run processing.
  - Activate live processing mode by simply adding the ``--live-processing`` argument and specifying the MinKNOW output directory as the input FAST5 directory.
  - Megalodon will continue to search for FAST5s until the ``final_summary*`` file is created by MinKNOW, indicating data production has completed.
- ``--outputs``

  - Specify desired outputs.
  - Options are ``basecalls``, ``mod_basecalls``, ``mappings``, ``variant_mappings``, ``mod_mappings``, ``per_read_variants``, ``per_read_mods``, ``variants``, and ``mods``.

    - ``mod_basecalls`` are output in a BAM file via the ``Mm`` and ``Ml`` tags `described by hts-specs here <https://github.com/samtools/hts-specs/pull/418>`_.
    - ``variant_mappings`` are intended for obtaining highly accurate phased variant genotypes, but also provide a nice genome browser visualiztion of per-read variant calls.

      - These mappings contain reference sequence at all positions except for per-read called variants. The base quality scores encode the likelihood for that reference anchored variant for use in the ``whathap`` phasing algorithm.
    - ``mod_mappings`` provide reference-anchored per-read modified base calls.

      - As of version 2.2, the default output uses the ``Mm`` and ``Ml`` hts-specs tags (see above) with all modified bases in one output file.
      - Specify the ``--mod-map-emulate-bisulfite`` option to output one BAM per modified base with modified bases converted using ``--mod-map-base-conv``

        - This file is useful for visualizing per-read modified base calls (e.g. IGV bisulfite mode for CpG calls).
        - This file may also allow a port to standard bisulfite pipelines that are capable of processing long-reads.
  - Default output is ``basecalls`` only.
- ``--output-directory``

  - Specify the directory to output results.
    Default ``megalodon_results``
- ``--overwrite``

  - Overwrite the ``--output-directory`` if it exists.
  - Note that this is a recursive file deletion and should be used with caution.

-----------------
Mapping Arguments
-----------------

- ``--mappings-format``

  - Format for ``mapping`` output.
  - Options include ``bam`` (default), ``cram``, and ``sam``.
  - As of version 2.2, mappings are no longer sorted by default.

    - Set ``--sort-mappings`` to sort mappings. If ``samtools`` is not in ``$PATH`` provide path to executable via the ``--samtools-executable`` argument.
- ``--reference``

  - Reference genome or transcriptome in FASTA or minimap2 index format.

    - If ``--reference`` is a minimap2 index and ``--mapping-format`` is ``cram``, provide FASTA reference via ``--cram-reference``.

--------------------------
Sequence Variant Arguments
--------------------------

- ``--haploid``

  - Compute sequence variants assuming a haploid reference. Default: diploid
- ``--variant-filename``

  - File containing putative variants in VCF/BCF format.

    - Variants file must be sorted.
    - If variant file is not compressed and indexed this will be performed before further processing.
  - Variants must be matched to the ``--reference`` provided.

-----------------------
Modified Base Arguments
-----------------------

- ``--mod-motif``

  - Restrict modified base results to the specified motifs.
  - This argument takes 3 values representing:

    1. Modified base single letter codes (see ``megalodon_extras modified_bases describe_alphabet`` command)
    2. Canonical sequence motif (may contain `ambiguity codes <https://droog.gs.washington.edu/parc/images/iupac.html>`_)
    3. Relative position (0-based) of the modified base within the canonical sequence motif
  - Multiple ``--mod-motif`` arguments can be provided to a single ``megalodon`` command.
  - If not provided (and ``per_read_mods`` or ``mods`` outputs requested) all relevant sites are tested (e.g. all ``C`` bases for ``5mC``).

    - Note that restricting to motifs of interest can save computationally expensive steps and is considered more than a simple post-processing filter.

--------------------------
Compute Resource Arguments
--------------------------

- ``--processes``

  - Number of CPU read-processing workers to spawn.
- ``--devices``

  - GPU devices to use for basecalling acceleration.
  - If not provided CPU basecalling will be performed.
  - Device names can be provided in the following formats: ``0``, ``cuda0`` or ``cuda:0``.
  - Multiple devices can be specified separated by a space.
