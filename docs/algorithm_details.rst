***************************
Megalodon Algorithm Details
***************************

This page describes the details of how megalodon processes the raw nanopore signal to produce highly-accurate modified base and SNP calls.

--------------------
Reference Aanchoring
--------------------

Megalodon's functionality centers on the anchoring of high-information neural network basecalling output to a reference sequence.
Given anchored neural network output, alternatives to the reference (either modified bases or canonical bases) are proposed and scored to produce the highest accuracy results.

The neural network output is anchored to the reference via standard read mapping of produced basecalls (maintaining the link to the neural network outputs) to the reference sequence.
If no reference mapping is produced (using ``minimap2`` via the ``mappy`` python interface) that read is not processed further (basecalls will be output if requested).
This standard read mapping is processed to produce a matching of each basecall with a reference position.
Positions within an insertion or deletion are left justified to the previous valid mapping position.
This constitutes the reference anchoring used for modified base and SNP calling steps.

-----------
SNP Calling
-----------

At each proposed SNP location the reference mapping identifies the range of the neural network output containing the sequence of interest.
The width of this sequence of interest is defined by the ``--snp-context-bases`` argument (specified individually for single base SNVs and indels; defaults 10 and 30 respectively).

The sequence scoring function performs the forward algorithm over the neural network output to produce a score for the reference and proposed alternative sequence.
The difference between these two scores is the assigned score for the proposed variant.
Lower score are evidence for the alternative sequence and positive scores are evidence for the reference sequence.

---------------------
Modified Base Calling
---------------------
