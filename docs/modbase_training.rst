**************************************
Megalodon Modified Base Model Training
**************************************

This page describes the process to generate modified base training data and train a basecalling model capable of also detecting modified bases.
The options provided from Megalodon for modified base training data annotation are each linked to a specific sample type.
If the provided options for do not perform sufficiently for a particular sample type, please open an issue on the `Megalodon issues page <https://github.com/nanoporetech/megalodon/issues>`_ with details for the sample type and intended training procedure.


------------------------------
Modified Base in Known Context
------------------------------

Given a sample with all locations matching a particular sequence motif modified (e.g. bacterial methylation) the ``--ref-mods-all-motifs`` option can be specified.
The ``--ref-mods-all-motifs`` argument takes ``4`` arguments.
These values are:

1. Single letter code for the modified base
2. Long name for the modified base
3. Sequence motif (using `IUPAC ambiguity codes <https://genome.ucsc.edu/goldenPath/help/iupac.html>`_)
4. Relative position for the modified base within the sequence motif (0-based coordinate)

The first two values are simply stored in the training data file.
When used for model training, these values will be saved for annotation of output files where appropriate (e.g. by Megalodon or Guppy output files).

Using this mode of modified base annotation, the specified basecalling model can be either a standard (canonical bases only) model or a modified base model.
Thus this method of annotation can be specified for modifications with no previous modeling (as long as the current basecalls map to the provided reference).
If the basecalling model specified is a modified base model, the single letter, long name, and corresponding canonical base attributes must be in agreement with the basecalling model.

As an example a native E. coli K12 sample contains 5mC methylation at a single motifs, ``CCWGG`` at the second position, and 6mA methylation at several motifs, ``GATC`` at the second position, ``AACNNNNNNGTGC`` at the second position, and ``GCACNNNNNNGTT`` at the third position (`source <https://www.neb.com/tools-and-resources/usage-guidelines/dam-and-dcm-methylases-of-e-coli>`_).
In order to annotate all of these modifications in a training data set the following command would be run:

::

   megalodon ecoli_k12_fast5s/ \
       --outputs signal_mappings \
       --reference ecoli_k12_reference.fa \
       --devices 0 --processes 40 \
       --ref-mods-all-motifs Z 5mC CCWGG 1 \
       --ref-mods-all-motifs Y 6mA GATC 1 \
       --ref-mods-all-motifs Y 6mA AACNNNNNNGTGC 1 \
       --ref-mods-all-motifs Y 6mA GCACNNNNNNGTT 2

Note that the single letter codes, ``Z`` and ``Y``, are arbitrary and can be set to any value desired by the user.
These values will be stored in the trained model for outputs where appropriate.

----------------------------
Modified Base Model Training
----------------------------

Given the above modified base annotated mapped signal file, a new modified base model can be trained with `Taiyaki <https://github.com/nanoporetech/taiyaki>`_.
Below is an example command to train a modified base model from the data prepared above.

::

   train_flipflop.py ./taiyaki/models/mLstm_cat_mod_flipflop.py \
       megalodon_results/signal_mappings.hdf5 --device 0
   # dump model to json format for use by guppy
   dump_json.py training/model_final.checkpoint \
       --output model_final.jsn


---------------------
Megalodon Calibration
---------------------

In order to produce the most accurate modified base calls, Megalodon requires the computation of a calibration file.
A ground truth sample containing known modified reference sites as well as known canonical base sites is required.
Similarly to sequence variant calibration, a modified base calibration file is created by running ``megalodon/scripts/generate_ground_truth_mod_llr_scores.py`` followed by ``megalodon/scripts/calibrate_mod_llr_scores.py``.

The ground truth modified base positions can be provided in one of two ways.

The first is by providing a single Megalodon results directory along with a CSV file containing ground truth modified and canonical sites.
The ``megalodon/scripts/create_mod_ground_truth.py`` script is provided to generate a ground truth CSV file from a set of bedmethyl files.
This first method is useful given a sample with variable methylation across the genome with alternative methylation evidence (e.g. human or plant methylation with bisulfite data).
If sites provided are strand specific (no ``--strand-offset`` specified to the ``create_mod_ground_truth.py`` script) specify the ``--strand-specific-sites`` option.

The second method to provide a modified base ground truth is via a control sample.
Using this method all called modified base sites in the main Megalodon results directory are assumed to be modified and all called sites in the ``--control-megalodon-results-dir`` directory are assumed to be canonical bases.
Called sites include only those specified via the ``--mod-motif`` argument to the main ``megalodon`` command.

Calibration files are computed from the database files, so there is no need to output per-read text files for calibration.

Using the generated ground truth statistics file, the ``calibrate_mod_llr_scores.py`` will then compute a calibration file for Megalodon use alongside this trained model.
It is highly recommended that the ``--out-pdf`` option is specified and inspected to ensure valid model results/performance.

Below an example of the calibration modified base plot is shown.
The top panel of this figure shows the distribution of modified base scores for the ground truth modified and canonical base per-read calls.
The distributions are modified slightly to force monotonicity to the peak of each distribution (original distributions are included in the plot).
These distributions should show reasonable separation for a well trained model.

The second panel shows the empirical (from ground truth data) and theoretical (from neural network scores).
Note that these distributions often differ quite significantly.
The theoretical/raw scores are what would be used given the ``--disable-mod-calibration`` option.
This underscores the importance of the calibration step.

The final panel shows the conversion from theoretical log likelihood ratio to emperical log likelihood ratio to be saved in the calibration file.
Note that monotonic smoothing is applied to this function as well.

The vertical bars on these plots, show the effect of common threshold choices on the ground truth data.

TODO: Add calibration figure


----------------------------------
Bootstrap Modified Base Annotation
----------------------------------

Once a modified base model is trained and calibration file computed, further models can be trained by marking data with this model.

.. warning::

  Great care should be taken when training a modified base basecalling model, especially with this method.
  The accuracy of refernce modified base markup is strongly indicative of the final modified base detection performance for a trained model.

In order to annotate modified bases using a modified base model the following command would be run:

::

   megalodon native_human_fast5s/ \
       --outputs signal_mappings \
       --reference human_reference.fa \
       --devices 0 --processes 40 \
       --mod-motif Z CG 0 \
       --mod-motif Y A 0 \
       --ref-include-mods \
       --guppy-params "-m model_final.jsn"

In this example 5mC calls are limited to CpG sites while 6mA calls are made at any motif.
Thus no 5mC sites will be marked in the output signal mapping file outside of the ``CG`` context.

The ``--ref-mod-threshold`` argument is provided to adjust the annotation based on modeling results.
By default the threshold to annotate a base as modified is a log likelihood ratio of ``0`` (i.e. the modified base is more likely than the canonical base based on emperically calibrated statistics).
In some samples this value may not be optimal.
The ``scripts/compute_mod_thresh_score.py`` command is provided for assistance in determining a reasonable value for this parameter.
