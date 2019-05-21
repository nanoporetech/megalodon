****************************
Advanced Megalodon Arguments
****************************

---------------
Model Arguments
---------------

- ``--flappie-model-name``

  -

----------------
Output Arguments
----------------

- ``--basecalls-format``

  - Currently only FASTA format is supported, but the option is provided for potential future capacity for alternative formats.

- ``--num-reads``

  - Number of reads to process. Intended for test runs on a subset of a dataset.

-----------------
Mapping Arguments
-----------------

- ``--prepend-chr-ref``

  - Prepend the ``chr`` string to reference chromosome names. Generally used for matching with VCF reference names.

-------------
SNP Arguments
-------------

- ``--disable-snp-calibration``

  -
- ``--heterozygous-factors``

  -
- ``--max-snp-size``

  -
- ``--prepend-chr-vcf``

  -
- ``--snp-all-paths``

  -
- ``--snp-calibration-filename``

  -
- ``--snp-context-bases``

  -
- ``--write-vcf-llr``

  -

-----------------------
Modified Base Arguments
-----------------------

- ``--disable-mod-calibration``

  -
- ``--mod-binary-threshold``

  -
- ``--mod-calibration-filename``

  -
- ``--mod-context-bases``

  -
- ``--mod-all-paths``

  -

----------------
Reference Output
----------------

This output category is intended for use in creating modified base training dataset for taiyaki training.

- ``--output-per-read-references``

  -
- ``--refs-include-mods``

  -
- ``--refs-include-snps``

  -
- ``--refs-percent-identity-threshold``

  -
- ``--refs-percent-coverage-threshold``

  -
- ``--refs-length-range``

-----------------------
Miscellaneous Arguments
-----------------------

- ``--database-safety``

  -
- ``--edge-buffer``

  -
- ``--not-recursive``

  -
- ``--suppress-progress``

  -
