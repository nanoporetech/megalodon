*****************************
``megalodon_extras validate``
*****************************

The ``megalodon_extras validate`` command group contains commands to validate mapping and modified base outputs from Megalodon.
Note that scripts to validate sequence variants are not provided here.
Other tools including `vcfeval <https://github.com/RealTimeGenomics/rtg-tools>`_ and `hapy.py <https://github.com/Illumina/hap.py>`_ are recommended for validation of sequence variant results.

-------------------------------------
``megalodon_extras validate results``
-------------------------------------

Validate per-read mapping and modified base results.

This command produces text and graphical summaries of mapping and modified base performance.

Mapping results include distributional statistics for each sample provided (output determined by ``--out-filename`` default ``stdout``), as well as a plot showing the distribution of mapping accuracies for each sample (see ``--out-pdf``).

----

.. figure::  _images/mapping_validate_results.png
   :align: center
   :width: 600

   Example ``validate results`` per-read mapping plot.

----

Per-read modified base results require a per-read ground truth for modified and canonical bases.
This can be provided by either 1) supplying a control sample via the ``--control-megalodon-results-dirs`` argument (assumes all modified base calls at ``--valid-sites`` in main Megalodon results are modified) or 2) providing a ground truth set of sites containing modified and canonical bases via the ``--ground-truth-data`` argument.
See the ``megalodon_extras modified_bases create_ground_truth`` command for help generating a ground truth file.

Per-read modified base results are analyzed to produce several metrics inlcuding the optimal `F1-score <https://en.wikipedia.org/wiki/F1_score>`_, `mean average precision <https://en.wikipedia.org/wiki/Evaluation_measures_(information_retrieval)#Mean_average_precision>`_ and `ROC AUC <https://en.wikipedia.org/wiki/Receiver_operating_characteristic>`_ amongst others.
By default, modified and canonical ground truth sites are filtered to contain the same number of statistics for these statistic computations.
It is highly recommended that this not be changed (via ``--allow-unbalance-classes``) as class imbalance can have a large effect on the statistics, thus effecting their comparison between runs and/or models.
Below are example graphical representations produced for per-read modified base validation.

----

.. figure::  _images/mod_pr_validate_results.png
   :align: center
   :width: 600

   Example ``validate results`` per-read modified base precision-recall curve plot.

.. figure::  _images/mod_roc_validate_results.png
   :align: center
   :width: 600

   Example ``validate results`` per-read modified base ROC curve plot.

.. figure::  _images/mod_dist_validate_results.png
   :align: center
   :width: 600

   Example ``validate results`` per-read modified base score distribution plot.

----

-------------------------------------------------------
``megalodon_extras validate aggregated_modified_bases``
-------------------------------------------------------


----------------------------------------------------
``megalodon_extras validate compare_modified_bases``
----------------------------------------------------
