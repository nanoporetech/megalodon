****************************
Advanced Megalodon Arguments
****************************

----------------
Output Arguments
----------------

- ``--basecalls-format``

  - Currently only FASTA format is supported, but the option is provided for potential future capacity for alternative formats.
- ``--num-reads``

  - Number of reads to process. Intended for test runs on a subset.
- ``--read-ids-filename``

  - A file containing ``read_ids`` to process (one per line).
  - Used in the variant phasing pipeline.

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
  - Two factors must be provided for single base subsititution variants and indels.
- ``--max-indel-size``

  - Maximum indel size to include in testing. Default: 50
- ``--variant-all-paths``

  - Compute the forward algorithm all paths score.
  - Default: Viterbi best-path score.
- ``--variant-calibration-filename``

  - File containing emperical calibration for sequence variant scores.
  - As created by megalodon/scripts/calibrate_variant_llr_scores.py.
  - Default: Load default calibration file.
- ``--variant-context-bases``

  - Context bases for single base SNP and indel calling. Default: [15, 30]
- ``--variant-locations-on-disk``

  - Force sequence variant locations to be stored only within on disk database table. This option will reduce the RAM memory requirement, but may drastically slow processing. Default: Store locations in memory and on disk.
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
- ``--mod-all-paths``

  - Compute forwards algorithm all paths score for modified base calls.
  - Default: Viterbi best-path score.
- ``--mod-aggregate-method``

  - Modified base aggregation method.
  - Choices: expectation_maximization (default), binary_threshold

- ``--mod-binary-threshold``

  - Hard threshold for modified base aggregation (probability of modified/canonical base).

    - Sites where no canonical or modified base acheives this level of confidence will be ignored in aggregation.
  - Default: 0.75
- ``--mod-calibration-filename``

  - File containing emperical calibration for modified base scores.
  - As created by megalodon/scripts/calibrate_mod_llr_scores.py.
  - Default: Load default calibration file.
- ``--mod-context-bases``

  - Context bases for modified base calling.
  - Default: 15

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
      - Modified base proportion estimates are stored in genotype fields specified by the single letter modified base encodings (definied in the model file).

- ``--mod-positions-on-disk``

  - Force modified base positions to be stored only within on disk database table. This option will reduce the RAM memory requirement, but may drastically slow processing. Default: Store positions in memory and on disk.
- ``--write-mod-log-probs``

  - Write per-read modified base log probabilities out in non-standard VCF field.

    - The ``LOG_PROBS`` field will contain semi-colon delimited log probabilities for modified base within each read at this site.
    - For sites with multiple modified bases, per-read calls for each modification type are separated by a comma as specified by the ``A`` genotype field type.

      - The order is consistent within each modification type so that per-read probabilities across all modification types can be reconstructed.

----------------
Reference Output
----------------

This output category is intended for use in generating reference sequences for taiyaki basecall model training.

- ``--output-per-read-references``

  - Flag to trigger this output type (similar to adding an option to ``--outputs``)
- ``--refs-include-mods``

  - Include modified base calls in per-read reference output.
- ``--refs-include-variants``

  - Include sequence variant calls in per-read reference output.
- ``--refs-percent-identity-threshold``

  - Only include reads with higher percent identity in per-read reference output.
- ``--refs-percent-coverage-threshold``

  -  Only include reads with higher read alignment coverage in per-read reference output.
- ``--refs-length-range``

  - Only include reads with specified read length in per-read reference output.

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
  - Default: 100
- ``--not-recursive``

  - Only search for fast5 read files directly found within the fast5 directory.
  - Default: search recursively
- ``--suppress-progress``

  - Suppress progress bar output.
