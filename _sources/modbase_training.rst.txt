**************************************
Megalodon Modified Base Model Training
**************************************

This page describes the process to generate modified base training data and train a basecalling model capable of also detecting modified bases.
The options provided from Megalodon for modified base training data annotation are each linked to a specific sample type.
If the provided options for do not perform sufficiently for a particular sample type, please open an issue on the `Megalodon issues page <https://github.com/nanoporetech/megalodon/issues>`_ with details for the sample type and intended training procedure.

Currently two sample types are supported for modified base training data markup:

1. Modified base in known sequence context (e.g. bacterial methylation)
2. Native (fractionally modified) sample with existing modified base basecalling model

Given the first sample above, a model will never be presented a modified and canonical base in close proximity (within the same training chunk).
Thus the second sample type is often required to produce a final model with sufficient performance across biologically relevant samples.

The highest accuracy markup of the native sample is imperative to achieve the highest performance modified base basecalling model.
In version 2.3, an adaptation was added to allow for markup informed by a reference methylation ground truth (fraction of modified reads at each reference site).
This ground truth can come from a number of sources and technologies.
Find more details below on how to use this method.

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
When used for model training, these values will be saved for annotation of output files where appropriate (e.g. Megalodon or Guppy output files).

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
       --ref-mods-all-motifs m 5mC CCWGG 1 \
       --ref-mods-all-motifs a 6mA GATC 1 \
       --ref-mods-all-motifs a 6mA AACNNNNNNGTGC 1 \
       --ref-mods-all-motifs a 6mA GCACNNNNNNGTT 2

Note that the single letter codes, ``m`` and ``a``, can be set to any value desired by the user.
It is recommended that the values follow `specifications found in hts-specs <https://github.com/jkbonfield/hts-specs/blob/methylation/SAMtags.tex#L528>`_.
These values will be stored in the trained model for outputs where appropriate.

----------------------------------
Bootstrap Modified Base Annotation
----------------------------------

Once a modified base model is trained (see above) and calibration file computed (see below), further models can be trained by marking modified base training data with this model.

.. warning::

  Great care should be taken when training a modified base basecalling model, especially with this method.
  The accuracy of reference modified base markup is strongly indicative of the final modified base detection performance for a trained model.

The following example assumes a trained model to detect 5mC (``m``) in ``CG`` sequence contexts (model specified in ``model_final.cfg``).
In order to annotate 5mC sites in a modified base training data set (``signal_mapping``) using a modified base model the following command would be run:

::

   megalodon native_human_fast5s/ \
       --outputs signal_mappings \
       --reference reference.fa.mmi \
       --devices 0 --processes 40 \
       --mod-motif m CG 0 \
       --ref-include-mods \
       --guppy-config model_final.cfg

The ``--ref-mod-threshold`` argument is provided to adjust the annotation based on modeling results.
By default the threshold to annotate a base as modified is a log likelihood ratio of ``0`` (i.e. the modified base is more likely than the canonical base based on empirically calibrated statistics).
In some samples this value may not be optimal.
The ``megalodon_extras modified_bases estimate_threshold`` command is provided for assistance in determining a reasonable value for this parameter.

-----------------------------------------------------
Ground Truth Aided Bootstrap Modified Base Annotation
-----------------------------------------------------

To further improve modified base training data markup, a reference anchored ground truth can be leveraged.
This method sets the modified base markup threshold at each reference position, informed by the provided ground truth.
This is similar to the ``--ref-mod-threshold`` argument, but this threshold is global to all reference positions.

The first step in this method is to call per-read modified bases on the native sample of interest.
This sample should contain sufficient depth, such that the identified modified base threshold at each position will be robust.
50X coverage is a rough target for sufficient coverage with this method.

Given a completed Megalodon run (``megalodon_results`` with ``--outputs mappings per_read_mods``) and the ground truth bedmethyl file (``ground_truth_methylation.CG.bed``) the following command will compute per-site modified base thresholds and low coverage sites.

::

   # Compute per-site thresholds
   megalodon_extras \
       modified_bases per_site_thresholds \
       megalodon_results \
       ground_truth_methylation.CG.bed \
       --strand-offset 1 \
       --ground-truth-coverage-pdf gt_cov.CG.pdf \
       --ground-truth-cov-min 50 \
       --nanopore-cov-min 50 \
       --out-blacklist-sites low_coverage_sites.CG.bed \
       --out-per-site-mod-thresholds site_mod_thresholds.CG.bed
   # sort low coverage sites for faster bedtools filtering
   sort -S 25% --parallel=56 -T /tmp/ \
       -k1,1V -k2,2n \
       -o low_coverage_sites.CG.sorted.bed
   # filter and sort first round mappings (% ident>90, read coverage>90%, length between 1000 and 20000)
   awk '$2 > 90 && $7 > 90 && $3 - $6 > 1000 && $3 - $6 < 20000 {print $8"\t"$10"\t"$11"\t"$1"\t.\t"$9}' \
       megalodon_results/mappings.summary.txt | \
       sort -S 25% --parallel=56 -T /tmp/ -k1,1V -k2,2n > \
       mappings.filtered.sorted.bed
   intersectBed \
       -a mappings.filtered.sorted.bed \
       -b low_coverage_sites.CG.sorted.bed -s -sorted -v | \
       awk '{print $4}' > train_read_ids.txt

Finally the ground truth aided modified base markup training data set is produced with the following command.

::

   megalodon \
       native_human_fast5s/ \
       --reference reference.fa.mmi \
       --output-directory per_site_markup_mega_res \
       --outputs per_read_mods signal_mappings \
       --mod-motif m CG 0 \
       --guppy-config model_final.cfg \
       --processes 40 --devices 0 \
       --ref-include-mods \
       --mod-per-site-threshold site_mod_thresholds.CG.bed \
       --read-ids-filename train_read_ids.txt

----------------------------
Modified Base Model Training
----------------------------

Given any of the above modified base annotated mapped signal files, a new modified base model can be trained with `Taiyaki <https://github.com/nanoporetech/taiyaki>`_.
Below is an example command to train a modified base model from the data prepared above and convert the final model for use with Guppy or Megalodon.

::

   train_flipflop.py ./taiyaki/models/mLstm_cat_mod_flipflop.py \
       megalodon_results/signal_mappings.hdf5 --device 0
   # dump model to json format for use by guppy
   dump_json.py training/model_final.checkpoint \
       --output model_final.jsn

The produced model should be referenced from a new Guppy config file.
The easiest way to obtain this would be to copy and edit the closest existing Guppy config file in the ``data`` directory of Guppy.

---------------------
Megalodon Calibration
---------------------

In order to produce the most accurate aggregated modified base calls, Megalodon requires the computation of a calibration file.
A ground truth sample containing known modified reference sites as well as known canonical base sites is required.
This can be the same as the model training data.
A modified base calibration file is created with the ``megalodon_extras calibrate generate_modified_base_stats`` and ``megalodon_extras calibrate modified_bases`` commands.
Please see the `calibration documentation page <https://nanoporetech.github.io/megalodon/extras_calibrate.html>`_ for more details about this process.
