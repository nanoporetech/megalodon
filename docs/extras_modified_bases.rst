***********************************
``megalodon_extras modified_bases``
***********************************

The ``megalodon_extras modified_bases`` command group contains various commands related to modified base processing within Megalodon.

------------------------------------------------------
``megalodon_extras modified_bases estimate_threshold``
------------------------------------------------------

Estimate the optimal global modified base score threshold such that the estimated proportion of bases modified is achieved when all sites are called.
This command is useful when producing the ``signal_mappings`` or ``per_read_refs`` outputs with modified bases annotated for basecaller training with Taiyaki.

This command works by estimating the proportion of bases modified from a sample using the most extreme calls (default most extreme 8%; set with ``--mod-percentile`` option) as a truth set.
This method assumes that the distribution of modified base and canonical scores is approximately balanced.
This may not be true for all models.
The plots produced by the ``megalodon_extras calibrate modified_bases`` command can assist in making this determination.

---------------------------------------------------
``megalodon_extras modified_bases update_database``
---------------------------------------------------

The Megalodon modified base database schema was updated at the release of Megalodon 2.0.
This command is provided in order to update a database from Megalodon <2.0 for use with Megalodon >=2.0.

-------------------------------------------------------
``megalodon_extras modified_bases create_ground_truth``
-------------------------------------------------------

Create a set of ground truth sites from a bedmethyl file.
This is a convenience command to apply a threshold to observed fractions of modified bases from bedmethyl files produced by Megalodon or other software.

The ``--strand-offset`` is provided to allow calls on opposite strands to be combined.
For example forward and reverse strand CpG calls can be merged by setting ``--strand-offset 1`` since the reverse strand position ``1`` downstream of the forward strand position correspond to the same biological methylation event.
