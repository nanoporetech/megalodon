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

- ``--outputs``

  - Specify desired outputs.
  - Options are ``basecalls``, ``mod_basecalls``, ``mappings``, ``variant_mappings``, ``mod_mappings``, ``per_read_variants``, ``per_read_mods``, ``variants``, and ``mods``.

    - ``mod_basecalls`` are currently output in an HDF5 file with a data set corresponding to each read (accessed via the ``read_id``).

      - This format is likely to change in the future.
    - ``variant_mappings`` are intended only for obtaining highly accurate phased variant genotypes.

      - These mappings contain reference sequence at all positions except for per-read called variants. The base quality scores encode the likelihood for that reference anchored variant for use in the whathap phasing algorithm.
      - This file is useful for visualizing per-read variant calls as well as potential variant phasing.
    - ``mod_mappings`` provide refernce-anchored per-read modified base calls.

      - These mappings contain the mappedreference sequence annotated with modified base calls at all instances of ``--mod-motif``s.
      - This file is useful for visualizing per-read modified base calls (e.g. IGV bisulfite mode for CpG calls).
      - This file may also allow a port to standard bisulfite pipelines that are capable of processing long-reads.
  - Default output is basecalls only.
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
  - Options include ``bam`` (defualt), ``cram``, and ``sam``.
- ``--reference``

  - Reference genome or transcriptiome in FASTA format.

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
- ``--variant-calibration-filename``

  - File containing emperical calibration for sequence variant scores.
  - As created by the ``megalodon_extras calibrate variants`` command.
  - Default: Load default calibration file for guppy config.

-----------------------
Modified Base Arguments
-----------------------

- ``--mod-motif``

  - Restrict modified base results to the specified motifs.
  - If not provided (and ``per_read_mods`` or ``mods`` outputs requested) all relevant sites are tested (e.g. all ``C`` bases for ``5mC``).
- ``--mod-calibration-filename``

  - File containing emperical calibration for modified base scores.
  - As created by ``megalodon_extras calibrate modified_bases`` command.
  - Default: Load default calibration file for guppy config.

-----------------------
Miscellaneous Arguments
-----------------------

- ``--processes``

  - Number of CPU read-processing workers to spawn.
- ``--devices``

  - GPU devices to use for basecalling acceleration.
  - If not provided CPU basecalling will be performed.
  - Device names can be provided in the following formats: ``0``, ``cuda0`` or ``cuda:0``.
  - Multiple devices can be specified separated by a space.
- ``--verbose-read-progress``

  - Output dynamic updates to potential issues during processing.
