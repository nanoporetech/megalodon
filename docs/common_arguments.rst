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

--------------
Model Argument
--------------

- ``--taiyaki-model-filename``

  - `taiyaki <https://github.com/nanoporetech/taiyaki>`_ basecalling model checkpoint file
  - In order to identify modified bases a model trained to identify those modifications must be provided.

    - Train a new modified base model using taiyaki.
  - By default the model included with megalodon is used.

    - This model is applicable to MinION or GridION R9.4.1 flowcells.
    - This model is equivalent to the "high accuracy" model for MinKNOW/guppy.
    - This model includes modified bases 5mC and 6mA in biological contexts: 5mC in human (CpG) and E. coli (CCWGG) contexts and 6mA in E. coli (GATC) context.

  - Guppy JSON-format models can be converted to taiyaki checkpoints/models with the ``taiyaki/bin/json_to_checkpoint.py`` script for use with megalodon.

----------------
Output Arguments
----------------

- ``--outputs``

  - Specify desired outputs.
  - Options are ``basecalls``, ``mod_basecalls``, ``mappings``, ``whatshap_mappings``, ``per_read_variants``, ``per_read_mods``, ``variants``, and ``mods``.

    - ``mod_basecalls`` are currently output in an HDF5 file with a data set corresponding to each read (accessed via the ``read_id``).
    - ``whatshap_mappings`` are intended only for obtaining highly accurate phased variant genotypes.

      - These mappings contain reference sequence at all positions except for per-read called variants. The base quality scores encode the likelihood for that reference anchored variant for use in the whathap phasing algorithm.
      - This file is useful for visualizing per-read variant calls as well as potential variant phasing.
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
- ``--write-variants-text``

  - Output per-read variants in text format.

    - Output includes columns: ``read_id``, ``chrm``, ``strand``, ``pos``, ``ref_log_prob``, ``alt_log_prob``, ``var_ref_seq``, ``var_alt_seq``, ``var_id``
    - Log probabilities are calibrated to match observed log-likelihood ratios from ground truth samples.

      - Reference log probabilities are included to make processing mutliple alternative allele sites easier to process.
    - Position is 0-based

-----------------------
Modified Base Arguments
-----------------------

- ``--mod-motif``

  - Restrict modified base results to the specified motifs.
  - If not provided (and ``per_read_mods`` or ``mods`` outputs requested) all relevant sites are tested (e.g. all ``C`` bases for ``5mC``).
- ``--write-mods-text``

  - Output per-read modified bases in text format.

    - Output includes columns: ``read_id``, ``chrm``, ``strand``, ``pos``, ``mod_log_probs``, ``can_log_prob``, ``mod_bases``, ``motif``
    - Log probabilities are calibrated to match observed log-likelihood ratios from ground truth samples.

      - Canonical log probabilities are included to make processing mutliple modification sites easier to process.

        - Note that the included model does not model such sites, but megalodon is capable of handling these sites (e.g. testing for 5mC and 5hmC simultaneously is supported given a basecalling model).
    - ``motif`` includes the searched motif (via ``--mod-motif``) as well as the relative modified base position within that motif (e.g. ``CG:0`` for provided ``--mod-motif Z CG 0``).
    - Position is 0-based

-----------------------
Taiyaki Chunk Arguments
-----------------------

- ``--chunk-size``

  - Size of individual chunks to run as input to neural network.
  - Smaller size will result in faster basecalling, but may reduce accuracy.
- ``--chunk-overlap``

  - Overlap between adjacent chunks fed to baescalling neural network.
  - Smaller size will result in faster basecalling, but may reduce accuracy.
- ``--devices``

  - GPU devices to use for basecalling acceleration.
  - If not provided CPU basecalling will be performed.
  - A separate GPU process will be spawned for each CPU worker process requested (spread evenly over specified ``--devices``).

    - Each GPU process must load the model parameters (~450MB for the default high-accuracy model).
    - Extra headroom for chunks must be allowed as well.
    - ``1`` process per 0.6 GB of GPU memory is a good default.
- ``--max-concurrent-chunks``

  - Maximum number of concurrent chunks to basecall at once.
  - Allows a global cap on GPU memory usage.
  - Changes to this parameter do not effect resulting basecalls.

-----------------------
Miscellaneous Arguments
-----------------------

- ``--processes``

  - Number of CPU worker processes to spawn.
- ``--verbose-read-progress``

  - Output dynamic updates to potential issues during processing.
