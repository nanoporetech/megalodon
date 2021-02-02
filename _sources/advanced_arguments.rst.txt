****************************
Advanced Megalodon Arguments
****************************

----------------------
Guppy Backend Argument
----------------------

- ``--do-not-use-guppy-server``

  - Use alternative basecalling backend
  - Alternatives are:

    - FAST5: Read called sequence and full posterior data from fast5 files.

      - This is the default when ``--do-not-use-guppy-server`` is set.
      - Note that this option requires ``--post_out`` be set when running Guppy and may increase the fast5 file size by 5-10X.
    - Taiyaki: Use the Taiyaki package basecalling interface

      - This requires a Taiyaki installation (potentially with GPU settings).
      - Trigger this mode by setting the ``--taiyaki-model-filename`` option.
      - This is much slower than Guppy and is generally intended for experimental models with either layers or architectures not supported by Guppy.
- ``--guppy-params``

  - Extra guppy server parameters.
  - Main purpose for optimal performance based on compute environment.
  - Quote parameters to avoid them being parsed by megalodon.
- ``--guppy-server-port``

  - Guppy server port.
  - Default: ``auto``
- ``--reads-per-guppy-batch``

  - Number of reads to send to Guppy per batch within each worker processes.
  - Default: ``50``
- ``--guppy-timeout``

  - Timeout to wait for guppy server to call a single read in seconds.
  - Default: ``5.0``
- ``--list-supported-guppy-configs``

  - List guppy configs with sequence variant and (if applicable) modified base support.

----------------
Output Arguments
----------------

- ``--basecalls-format``

  - Select either ``fastq`` (default) or ``fasta`` format for basecalls output.
- ``--num-reads``

  - Number of reads to process. Intended for test runs on a subset.
- ``--read-ids-filename``

  - A file containing ``read_ids`` to process (one per line).
  - Used in the variant phasing pipeline.
- ``--mod-min-prob``

  - Only include modified base probabilities greater than this value in ``mod_basecalls`` and ``mod_mappings`` outputs.
  - Default: ``0.01`` (``1%``)

-----------------
Mapping Arguments
-----------------

- ``--cram-reference``

  - If ``--reference`` is a minimap2 index, the associated FASTA reference needs to be provided for ``--mappings-format cram``.
- ``--samtools-executable``

  - Samtools executable or path for sorting and indexing all mappings.
  - Default: ``samtools``
- ``--sort-mappings``

  - Perform sorting and indexing of mapping output files.
  - This can take considerable time for larger runs and thus is off by default.

--------------------------
Sequence Variant Arguments
--------------------------

- ``--context-min-alt-prob``

  - Minimum per-read variant probability to include a variant in second round of variant evaluation (including context variants).

- ``--disable-variant-calibration``

  - Use raw neural network sequence variant scores.
  - This option should be set when calibrating a new model.
  - Default: Calibrate scores as described in ``--variant-calibration-filename``
- ``--heterozygous-factors``

  - Bayes factor used when computing heterozygous probabilities in diploid variant calling mode.
  - Two factors must be provided for single base substitution variants and indels.
- ``--max-indel-size``

  - Maximum indel size to include in testing. Default: 50
- ``--variant-all-paths``

  - Compute the forward algorithm all paths score.
  - Default: Viterbi best-path score.
- ``--variants-are-atomized``

  - Input variants have been atomized (with ``megalodon_extras variants atomize``).
  - This saves compute time, but has unpredictable behavior if variants are not atomized.
- ``--variant-calibration-filename``

  - File containing empirical calibration for sequence variant scores.
  - As created by the ``megalodon_extras calibrate variants`` command.
  - Default: Load default calibration file for guppy config.
- ``--variant-context-bases``

  - Context bases for single base SNP and indel calling. Default: [15, 30]
- ``--variant-locations-on-disk``

  - Force sequence variant locations to be stored only within on disk database table. This option will reduce the RAM memory requirement, but may drastically slow processing. Default: Store locations in memory and on disk.
- ``--write-variants-text``

  - Output per-read variants in text format.

    - Output includes columns: ``read_id``, ``chrm``, ``strand``, ``pos``, ``ref_log_prob``, ``alt_log_prob``, ``var_ref_seq``, ``var_alt_seq``, ``var_id``
    - Log probabilities are calibrated to match observed log-likelihood ratios from ground truth samples.

      - Reference log probabilities are included to make processing multiple alternative allele sites easier to process.
    - Position is 0-based
- ``--write-vcf-log-probs``

  - Write per-read alt log probabilities out in non-standard VCF field.

    - The ``LOG_PROBS`` field will contain semi-colon delimited log probabilities for each read at this site.
    - For sites with multiple alternative alleles, per-read calls for each allele are separated by a comma as specified by the ``A`` genotype field type.

      - The order is consistent within each allele so that per-read probabilities across all alleles can be reconstructed.

-----------------------
Modified Base Arguments
-----------------------

- ``--disable-mod-calibration``

  - Use raw modified base scores from the network.
  - This option should be set when calibrating a new model.
  - Default: Calibrate scores as described in ``--mod-calibration-filename``
- ``--mod-aggregate-method``

  - Modified base aggregation method.
  - Choices: expectation_maximization (default), binary_threshold
- ``--mod-all-paths``

  - Compute forwards algorithm all paths score for modified base calls.
  - Default: Viterbi best-path score.
- ``--mod-binary-threshold``

  - Hard threshold for modified base aggregation (probability of modified/canonical base).

    - Sites where no canonical or modified base achieves this level of confidence will be ignored in aggregation.
  - Default: 0.75
- ``--mod-calibration-filename``

  - File containing empirical calibration for modified base scores.
  - As created by ``megalodon_extras calibrate modified_bases`` command.
  - Default: Load default calibration file for guppy config.
- ``--mod-database-timeout``

  - Timeout in seconds for modified base database operations.
  - Default: 5 seconds
- ``--mod-context-bases``

  - Context bases for modified base calling.
  - Default: 15
- ``--mod-map-emulate-bisulfite``

  - For ``mod_mappings`` output, emulate bisulfite output by converting called modified bases using "--mod-map-base-conv" argument.
  - As of version 2.2, the default ``mod_mappings`` output uses the ``Mm`` and ``Ml`` hts-specs tags (see above) with all modified bases in one output file.
- ``--mod-map-base-conv``

  - For ``mod_mappings`` output, convert called bases.

    - For example, to mimic bisulfite output use: ``--mod-map-base-conv C T --mod-map-base-conv Z C``
    - This is option useful since the BAM format does support modified bases and will convert all alternative bases to ``N``s for storage in BAM/CRAM format.
  - Note additional formats may be supported in the future once finalized in hts-specs.
- ``--mod-output-formats``

  - Modified base aggregated output format(s).
  - Default: ``bedmethyl``
  - Options: ``bedmethyl``, ``modvcf``, ``wiggle``

    - ``bedmethyl`` format produces one file per modification type.

      - This format is specified by the `ENCODE consortium <https://www.encodeproject.org/data-standards/wgbs/>`_.
    - ``modvcf`` is a slight variant to the VCF format used for sequence variant reporting.

      - This format produces a single file containing all modifications.
      - The format adds a ``SN`` info field as modified bases occur in a stranded manner unlike sequence variants (e.g. hemi-methylation).
      - A genotype field ``VALID_DP`` indicates the number of reads included in the proportion modified calculation.
      - Modified base proportion estimates are stored in genotype fields specified by the single letter modified base encodings (defined in the model file).
- ``--write-mod-log-probs``

  - Write per-read modified base log probabilities out in non-standard VCF field.

    - The ``LOG_PROBS`` field will contain semi-colon delimited log probabilities for modified base within each read at this site.
    - For sites with multiple modified bases, per-read calls for each modification type are separated by a comma as specified by the ``A`` genotype field type.

      - The order is consistent within each modification type so that per-read probabilities across all modification types can be reconstructed.
- ``--write-mods-text``

  - Output per-read modified bases in text format.

    - Output includes columns: ``read_id``, ``chrm``, ``strand``, ``pos``, ``mod_log_probs``, ``can_log_prob``, ``mod_bases``, ``motif``
    - Log probabilities are calibrated to match observed log-likelihood ratios from ground truth samples.

      - Canonical log probabilities are included to make processing multiple modification sites easier to process.

        - Megalodon is capable of handling multiple modified bases per site with appropriate model (e.g. testing for 5mC and 5hmC simultaneously is supported given a basecalling model).
    - ``motif`` includes the searched motif (via ``--mod-motif``) as well as the relative modified base position within that motif (e.g. ``CG:0`` for provided ``--mod-motif Z CG 0``).
    - Position is 0-based

-------------------------
Taiyaki Backend Arguments
-------------------------

- ``--chunk-size``

  - Size of individual chunks to run as input to neural network.
  - Smaller size will result in faster basecalling, but may reduce accuracy.
- ``--chunk-overlap``

  - Overlap between adjacent chunks fed to basecalling neural network.
  - Smaller size will result in faster basecalling, but may reduce accuracy.
- ``--max-concurrent-chunks``

  - Maximum number of concurrent chunks to basecall at once.
  - Allows a global cap on GPU memory usage.
  - Changes to this parameter do not effect resulting basecalls.
- ``--taiyaki-model-filename``

  - `taiyaki <https://github.com/nanoporetech/taiyaki>`_ basecalling model checkpoint file
  - In order to identify modified bases a model trained to identify those modifications must be provided.

    - Train a new modified base model using taiyaki.

  - Guppy JSON-format models can be converted to taiyaki checkpoints/models with the ``taiyaki/bin/json_to_checkpoint.py`` script for use with megalodon.

-------------------------------
Reference/Signal Mapping Output
-------------------------------

This output category is intended for use in generating reference sequences or signal mapping files for taiyaki basecall model training.

- ``--ref-include-mods``

  - Include modified base calls in ``per_read_refs`` or ``signal_mappings`` outputs.
- ``--ref-include-variants``

  - Include sequence variant calls in per-read reference output.
- ``--ref-length-range``

  - Only include reads with specified read length in per-read reference output.
- ``--ref-percent-identity-threshold``

  - Only include reads with higher percent identity in per-read reference output.
- ``--ref-percent-coverage-threshold``

  -  Only include reads with higher read alignment coverage in per-read reference output.
- ``--ref-mods-all-motifs``

  - Annotate all ``--mod-motif`` occurrences as modified.
  - Requires that `--ref-include-mods`` is set.
- ``--ref-mod-threshold``

  - Threshold (in ``log(can_prob/mod_prob)`` space) used to annotate a modified bases in ``signal_mappings`` or ``per_read_refs`` outputs.
  - See ``megalodon_extras modified_bases estimate_threshold`` command for help computing this threshold.
  - Requires that `--ref-include-mods`` is set.

--------------------------
Compute Resource Arguments
--------------------------

- ``--num-read-enumeration-threads``

  - Number of parallel threads to use for read enumeration.

    - This number of threads will be opened in a single read enumeration process and each signal extraction process (see next argument).
  - This value can be increased if the input queue remains empty.
  - Default: ``8``
- ``--num-extract-signal-processes``

  - Number of parallel processes to use for signal extraction.

    - Accessing data and metadata from FAST5 files requires some compute resources. For this reason, multiple processes must be spawned to achieve the highest performance on some systems.
  - This value can be increased if the input queue remains empty.
  - Default: ``2``

-----------------------
Miscellaneous Arguments
-----------------------

- ``--database-safety``

  - Setting for database performance versus corruption protection.

    - Options:

      - 0 (DB corruption on application crash)
      - 1 (Default; DB corruption on system crash)
      - 2 (DB safe mode)
- ``--edge-buffer``

  - Do not process sequence variant or modified base calls near edge of read mapping.
  - Default: 30
- ``--not-recursive``

  - Only search for fast5 read files directly found within the fast5 directory.
  - Default: search recursively
- ``--suppress-progress``

  - Suppress progress bar output.
- ``--suppress-queues-status``

  - Suppress dynamic status of output queues.
  - These queues are helpful for diagnosing I/O issues.
- ``--verbose-read-progress``

  - Output dynamic updates to potential issues during processing.
  - Default: ``3``
