************************
Megalodon Model Training
************************

This page describes how to use megalodon to prepare training data and train a new basecalling model using taiyaki.

.. note::

   Preparation of training data via megalodon requires a basecalling model that can produce valid reference mappings.
   If reference valid mappings cannot be produced for a set of reads, model training will not proceed successfully.

----------------
Data Preparation
----------------

Data preparation involves adding the ``--outputs signal_mappings`` argument to a megalodon call.
This will produce a ``signal_mappings.hdf5`` file in the specified megalodon output directory.
For each read producing a valid reference mapping, this file contains a mapping between the raw signal and the mapped reference bases.
This file can then be directly passed to the taiyaki ``train_flipflop.py`` command for model training.

::

   # set unix environment variables for optimal processing
   export OPENBLAS_NUM_THREADS=1
   export OMP_NUM_THREADS=1
   # run megalodon outputting signal mappings and read mappings
   megalodon raw_fast5s/ \
       --outputs mappings signal_mappings \
       --reference reference.fa \
       --devices 0 --processes 40

   # set environment variables for optimal taiyaki training
   export OMP_NUM_THREADS=40
   # run taiyaki training
   train_flipflop.py ./taiyaki/models/mLstm_flipflop.py \
       megalodon_results/signal_mappings.hdf5 --device 0

Once training completes, the ``training/model_final.checkpoint`` contains the model.
This can be converted to a guppy compatible model with the ``taiyaki/bin/dump_json.py`` script.
A guppy config with appropriate settings should also be produced for new models.

----------------------
Signal Mapping Options
----------------------

Several options are available to control the behavior of the ``signal_mappings`` output.

- ``--ref-length-range``

  - Only allow reads with a reference mapping length within this range into the output.
- ``--ref-percent-identity-threshold``

  - Only include reads with higher mappig percent identity in signal_mappings output.
- ``--ref-percent-coverage-threshold``

  - Only include reads with higher read alignment coverage in signal_mappings output.
- ``--ref-include-variants``

  - This option replaces the reference sequence with more likley proposed alternative sequences as called in the ``per_read_variants`` output.
  - Cannot specify both this option and ``--ref-include-mods``.

------------------------------
Modified Base Data Preparation
------------------------------

Several options are available to allow the output of modified bases into the output training data.
In order to output modified base training data a modbase model must be provided.

.. warning::

  Great care should be taken when training a modified base basecalling model.
  The accuracy of refernce modified base markup is strongly indicative of the final modified base detection performance for a trained model.

- ``--ref-include-mods``

  - Modified base calls will be included in signal mapping output.
  - If provided, calls will not be made outside of sequences specified by ``--mod-motif``.

    - Note that multiple ``--mod-motifs`` values may be specified.
    - By default, modified base calls will be made in all sequence contexts.
  - By default, the most likely modified base or canonical call will be included in the output reference where applicable.
  - If a modbase model is provided, but this option is not specified, the produced dataset will included the modified bases within the alphabet, but no modified bases will be annotated in the read training sequences.
- ``--ref-mods-all-motifs``

  - This option will mark all tested sites as modified regardless of the score at each tested reference site.
  - This option is applicable to samples where all instances of a particular motif are known to be modified (e.g. bacterial methylation).
- ``--ref-mod-threshold``

  - This option allows the user to manually set the log likelihood ratio threshold value for modidifed base reference markup.
  - See the ``scripts/compute_mod_thresh_score.py`` command for assistance on determining a reasonable value for this parameter.

---------------------
Megalodon Calibration
---------------------

TODO: Describe megalodon calibration.
