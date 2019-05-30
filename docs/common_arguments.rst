****************
Common Arguments
****************

------------------
Required Arguments
------------------

- ``fast5s_dir``

  - Path to directory containing raw FAST5-format nanopore reads.
  - Both single and mulit FAST5 formats are supported.
  - Default searches recursively for fast5 read files. To search only one-level specify `--not-recursive`.
- ``--taiyaki-model-filename``

  - `taiyaki <https://github.com/nanoporetech/taiyaki>`_ basecalling model checkpoint file

----------------
Output Arguments
----------------

- ``--outputs``

  -
- ``--output-directory``

  -
- ``--overwrite``

  -

-----------------
Mapping Arguments
-----------------

- ``--mappings-format``

  -
- ``--reference``

  -

-------------
SNP Arguments
-------------

- ``--haploid``

  -
- ``--variant-filename``

  -
- ``--write-snps-text``

  -

-----------------------
Modified Base Arguments
-----------------------

- ``--mod-motif``

  -
- ``--write-mods-text``

  -

-----------------------
Taiyaki Chunk Arguments
-----------------------

- ``--chunk_size``

  -
- ``--chunk_overlap``

  -
- ``--devices``

  -
- ``--max_concurrent_chunks``

  -

-----------------------
Miscellaneous Arguments
-----------------------

- ``--processes``

  -
- ``--verbose-read-progress``

  -
