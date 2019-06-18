****************
Common Arguments
****************

------------------
Required Arguments
------------------

- ``fast5s_dir``

  - Path to directory containing raw FAST5-format nanopore reads.
  - Both single and multi FAST5 formats are supported.
  - Default searches recursively for fast5 read files. To search only one-level specify `--not-recursive`.
- ``--taiyaki-model-filename``

  - `taiyaki <https://github.com/nanoporetech/taiyaki>`_ basecalling model checkpoint file
  - In order to identify modified bases a model trained to identify those modifications must be provided.

    - Train a new modified base model using taiyaki.

----------------
Output Arguments
----------------

- ``--outputs``

  - Specify desired outputs.
  - Options are ``basecalls``, ``mod_basecalls``, ``mappings``, ``per_read_snps``, ``per_read_mods``, ``snp``, and ``mods``.
  - Default is to output basecalls only.
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

-------------
SNP Arguments
-------------

- ``--haploid``

  - Compute sequence variants assuming a haploid reference. Default assume diploid.
- ``--variant-filename``

  - File containing putative variants in VCF/BCF format.
  - If variant file is not sorted, compressed and indexed this will be performed before further processing.
  - Variants must be matched to the ``--reference`` provided.
- ``--write-snps-text``

  - Output per-read SNPs in text format.

    - Output includes columns: ``read_id``, ``chrm``, ``strand``, ``pos``, ``scores``, ``snp_ref_seq``, ``snp_alt_seqs``, ``snp_id``
    - Scores are log probabilities of the alternative allele.

      - Probabilities are calibrated to match observed log-likelihood ratios from ground truth samples.
      - Note that for multi-allelic sites all alternative alleles must be extracted in order to obtain the reference probability.
    - Position is 0-based

-----------------------
Modified Base Arguments
-----------------------

- ``--mod-motif``

  - Restrict modified base results to the specified motifs.
  - If not provided (and ``per_read_mods`` or ``mods`` outputs requested) all relevant sites are tested.
- ``--write-mods-text``

  - Output per-read modified bases in text format.

    - Output includes columns: ``read_id``, ``chrm``, ``strand``, ``pos``, ``mod_log_probs``, ``mod_bases``, ``motif``
    - Log probabilities are comma-separated corresponding to modified bases in the next field

      - Log-probabilities are calibrated to match observed log likelihood ratios.
    - ``motif`` includes the reference sequence matching the specified motif (via ``--mod-motif``) as well as the relative modified base position within that motif
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
  - In not provided CPU basecalling will be performed.
  - A separate GPU process will be spawned for each CPU worker process requested (spread evenly over specified ``--devices``)..

    - Each GPU process must load the model parameters ~450MB for the default high-accuracy model).
    - Extra headroom for chunks must be allowed as well.
    - ``1`` process per GB of GPU memory is a good default.
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
