***********************************
``megalodon_extras modified_bases``
***********************************

The ``megalodon_extras modified_bases`` command group contains various commands related to modified base processing within Megalodon.

-----------------------------------------------------
``megalodon_extras modified_bases describe_alphabet``
-----------------------------------------------------

Describe the alphabet, including modified bases, found in a given model.

This command is useful to determine the syntax for specifying arguments related to modified base detection.

Note that originally modified bases were specified by arbitrary values, but recent and future models will attempt to follow single letter codes specified by `samtools hts-spec (currently an open issue) <https://github.com/samtools/hts-specs/pull/418>`_.

A minimal subset of the model specifications from the main ``megalodon`` command are available to specify the model exactly as in the main command.
The model will be loaded as normal and the alphabet used will be printed.

------------------------------------------------------
``megalodon_extras modified_bases estimate_threshold``
------------------------------------------------------

Estimate the optimal global modified base score threshold such that the estimated proportion of bases modified is achieved when all sites are called.
This command is useful when producing the ``signal_mappings`` or ``per_read_refs`` outputs with modified bases annotated for basecaller training with Taiyaki.

This command works by estimating the proportion of bases modified from a sample using the most extreme calls (default most extreme 8%; set with ``--mod-percentile`` option) as a truth set.
This method assumes that the distribution of modified base and canonical scores is approximately balanced.
This may not be true for all models.
The plots produced by the ``megalodon_extras calibrate modified_bases`` command can assist in making this determination.

----------------------------------------------------------------
[DEPRECATED] ``megalodon_extras modified_bases update_database``
----------------------------------------------------------------

This method has not been updated for the most recent modified base schema and thus is currently deprecated.

-------------------------------------------------------
``megalodon_extras modified_bases create_ground_truth``
-------------------------------------------------------

Create a set of ground truth sites from a bedmethyl file.
This is a convenience command to apply a threshold to observed fractions of modified bases from bedmethyl files produced by Megalodon or other software.

The ``--strand-offset`` is provided to allow calls on opposite strands to be combined.
For example forward and reverse strand CpG calls can be merged by setting ``--strand-offset 1`` since the reverse strand position ``1`` downstream of the forward strand position correspond to the same biological methylation event.

----------------------------------------------------
``megalodon_extras modified_bases create_motif_bed``
----------------------------------------------------

This is a helper command to take a reference FASTA file and produce a BED file with all the locations of a motif of interest.
This can be useful for a number of modified base pipelines.

-------------------------------------------------------
``megalodon_extras modified_bases per_site_thresholds``
-------------------------------------------------------

See the `ground truth aided bootstrap modified base annotation <https://nanoporetech.github.io/megalodon/modbase_training.html#ground-truth-aided-bootstrap-modified-base-annotation>`_ tutorial for more complete instructions on using this command.

This command is targeted at creating higher accuracy modified base training data sets (mapped signal files) from a preliminary modified base model and a fractionally (at each reference site) modified sample.
This command takes as input a Megalodon run with ``--outputs per_read_mods mappings`` and a ground truth bedmethyl file and produces a custom modified base threshold at each reference site.
These modified base thresholds will match the Megalodon modified base statistics with the ground truth faction of modified bases at each reference site.
The derived thresholds should then be supplied to Megalodon via the ``--mod-per-site-threshold`` argument along with the ``--outputs signal_mappings``.

This command will also output a low coverage BED file containing reference sites covered by either too few ground truth or nanopore reads.
This should be used to filter reads for the final training mapped signal data set.
A read covering any low coverage sites will not be marked up as accurately and thus should not be included in model training.

--------------------------------------------------
``megalodon_extras modified_bases index_database``
--------------------------------------------------

In certain instances a megalodon run may end unexpectedly (e.g. out of memory error).
In most cases the modified bases database is not corrupted by such an unexpected run termination.
This will leave the modified base database without having completed the indexing step which is required for most downstream uses.
This command will produce the index as a separate step in such instances.

--------------------------------------------------
``megalodon_extras modified_bases split_by_motif``
--------------------------------------------------

Split an input modified base database into smaller databases based on reference sequence motifs.
This command enables the computationally expensive ``megalodon`` command to be run only once, while still allowing motif specific analyses.
