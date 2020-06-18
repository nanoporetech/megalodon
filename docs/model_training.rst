************************
Megalodon Model Training
************************

This page describes how to use Megalodon to prepare training data and train a new basecalling model using Taiyaki.
For modified base data preparation and model training documentation see the :doc:`modified base training documentation </modbase_training>` page.

.. note::

   Preparation of training data via Megalodon requires a basecalling model that can produce valid reference mappings.
   If valid reference mappings using ``minimap2`` cannot be produced for a set of reads, model training will not proceed successfully.

----------------
Data Preparation
----------------

To produce a training data ("mapped signal") file the ``--outputs signal_mappings`` argument should be added to a Megalodon call.
This will produce a ``signal_mappings.hdf5`` file in the specified Megalodon output directory.
For each read producing a valid reference mapping, this file contains a mapping between the raw signal and the mapped reference bases.
This file can then be directly passed to the Taiyaki ``train_flipflop.py`` command for model training.

::

   # run megalodon; output signal mappings
   megalodon raw_fast5s/ \
       --outputs signal_mappings \
       --reference reference.fa \
       --devices 0 --processes 40

   # run taiyaki training
   train_flipflop.py ./taiyaki/models/mLstm_flipflop.py \
       megalodon_results/signal_mappings.hdf5 --device 0

Once training completes, the ``training/model_final.checkpoint`` contains the model.
This can be converted to a guppy compatible model with the ``taiyaki/bin/dump_json.py`` script.
A guppy config with appropriate settings should also be produced for new models.

.. note::

   For optimal performance, it is recommended that the ``OMP_NUM_THREADS`` unix environment variable be set to ``1`` for the above Megalodon command and a larger value for the Taiyaki training command.


----------------------
Signal Mapping Options
----------------------

Several options are available to control the behavior of the ``signal_mappings`` output.

- ``--ref-length-range``

  - Only allow reads with a reference mapping length within this range into the output.
- ``--ref-percent-identity-threshold``

  - Only include reads with higher mapping percent identity in signal_mappings output.
- ``--ref-percent-coverage-threshold``

  - Only include reads with higher read alignment coverage in signal_mappings output.
- ``--ref-include-variants``

  - This option replaces the reference sequence with more likely proposed alternative sequences as called in the ``per_read_variants`` output.
  - Cannot specify both this option and ``--ref-include-mods``.


---------------------
Megalodon Calibration
---------------------

When a new model is trained, the produced scores must be calibrated to achieve optimal aggregated results (over reads).
Once produced, calibration files can be passed to Megalodon via the ``--variant-calibration-filename`` and ``--mod-calibration-filename`` arguments.

Sequence variant calibration requires a ground truth against which to compute scores.
For sequence variants, a high quality reference for a set of reads will suffice for this requirement.
Random sequence variants are proposed and scored in order to create distributions over which to calibrate the produced scores.
In order to create a sequence variant calibration file, run ``megalodon/scripts/generate_ground_truth_variant_llr_scores.py`` followed by ``megalodon/scripts/calibrate_variant_llr_scores.py``.
The optional ``--out-pdf`` provides visualization of the likelihood ratio score correction.
