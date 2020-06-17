******************************
``megalodon_extras calibrate``
******************************

The ``megalodon_extras calibrate`` command group contains commands to produce Megalodon modified base and sequence variant calibration files for basecalling models.
When a new basecalling model is trained a calibration file must be produced in order to obtain the most accurate modified base and sequence variant calls.
Without a calibration file the ``--disable-mod-calibration`` or ``--disable-variant-calibration`` flags may be set, but aggregated results will likely be much less accurate.

Calibration file estimation is broken down into two step, 1) ground truth statistic generation (``megalodon_extras calibrate generate_*`` commands and 2) calibration estimation.
Detail of these steps for both modified base and sequence variants can be found below.

Note that the plots produced by the calibration procedure (with examples shown below) are stored in the repository for each released model in the GitHub repository (``megalodon/model_data/``).

---------------------------------------------
``megalodon_extras calibrate modified_bases``
---------------------------------------------

Estimate modified base calibration file.

Given a set of ground truth modified bases statistics, via ``--ground-truth-llrs`` argument, compute empirical probabilities of a modified base.
This command computes the empirical log-likelihood ratio over windows of observed modified base scores.
This process involves several steps to ensure certain characteristics of the generating distributions.
A separate calibration will be computed and stored in the output calibration file for each modified base found in the ground truth file.

These steps are visualized in the example plot below, which can be produced for any new calibration file by providing the ``--out-pdf`` argument.
The top facet of this plot shows the distribution of theoretical modified base log-likelihood ratios produced by the basecalling model.
These distributions are smoothed such that they are monotonic from either extreme to the peak of the densities.
The middle facet shows the inferred empirical probability that a base is modified given the theoretical modified base score produced by the basecaller.
The final facet shows the same probabilities, but in log-likelihood space.
A constraint is enforced on this function such that the value is monotonically increasing (red: before monotonic constraint; yellow: after monotonic constraint).
The three vertical lines indicate common threshold values for modified base aggregation.
Note that the fraction of data ignored at each threshold level is annotated in the figure legend.

----

.. figure::  _images/modified_base_calibration.png
   :align: center
   :width: 600

   Visualization of modified base calibration method.

----

---------------------------------------
``megalodon_extras calibrate variants``
---------------------------------------

Estimate sequence variant calibration file.

Given a set of ground truth sequence variant statistics, via ``--ground-truth-llrs`` argument, compute empirical probabilities of a sequence variant.
This command computes the empirical log-likelihood ratio over windows of observed sequence variant scores.
This process involves several steps to ensure certain characteristics of the generating distributions.
This procedure is largely the same as the modified base calibration step, but the variants are grouped into categories based on the type of ground truth sequence variant.
Note that the vertical bars are not present in these plots as sequence variant per-read statistics are combined in a probabilistic fashion and not based on a hard threshold.

----

.. figure::  _images/sequence_variant_calibration.png
   :align: center
   :width: 600

   Visualization of sequence variant calibration method.

----

-----------------------------------------------------------
``megalodon_extras calibrate generate_modified_base_stats``
-----------------------------------------------------------

Generate ground truth modified base statistics.

Extract modified base scores from a completed Megalodon run.
The ground truth can be provided by either 1) passing a control results directory (this assumes that the first Megalodon results are completely modified) or 2) provide a ``--ground-truth-data`` file containing modified and canonical base locations within a sample.
See the ``megalodon_extras modified_bases create_ground_truth`` command for help producing a ground truth CSV file.

-----------------------------------------------------
``megalodon_extras calibrate generate_variant_stats``
-----------------------------------------------------

Generate ground truth sequence variant statistics.

This method produces ground truth sequence variant statistics by proposing alternatives to a reference sequence.
It is thus assumed that the mapping location for each read contains the correct reference sequence.
It is advised to select a set of reads with high quality mappings to a high quality reference for the sample.

This command performs basecalling and read mapping as in the main Megalodon command.
Variants are then randomly proposed and scored for a random set of sites across each read.
"Correct" variants are not produced by default due to the computational overhead required to map full reads to the "incorrect" reference.
This functionality is provided on an experimental basis via the ``--compute-false-reference-scores`` flag, but these scores are not currently accepted by the ``megalodon_extras calibrate variants`` command.

---------------------------------------------------
``megalodon_extras calibrate merge_modified_bases``
---------------------------------------------------

Merge modified base calibration files.

In some cases the ground truth source for one modified base my come from a different source than another modified base.
In this case calibration files can be computed separately and combined with this command.
If multiple calibration files contain calibration for the same modified base, the calibration from the file listed first will be stored.
