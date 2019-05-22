***************************
Megalodon Algorithm Details
***************************

This page describes the details of how megalodon processes the raw nanopore signal to produce highly-accurate modified base and SNP calls.

------------
Base Calling
------------

Basecalling is performed as in taiyaki and guppy.
This involves taking raw nanopore signal, normalizing (using median and MAD scaling), chunking, running a recurrent neural network and decoding using the forward-backward algorithm.
These steps are described in the taiyaki documentation.

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

Megalodon currently filters alleles over a certain maximum size (default 5) as performance on larger indels has not currenty been tested.
An absolute maximum size of 27 is applied to to a current software reqstriction, but a fix has been designed.
Additionally, multi-allelic loci are currently split and treated as multiple bi-allelic variants.
Thus a heterozygous call including both alternative alleles is not currently capable from the this version of megalodon, but agian a fix has been designed for this issue.

At each proposed variant a region of context sequence around the variant is extracted.
The context sequence allows the scoring algorithm to traverse slightly different paths through the local neural network output.
The width of this sequence of interest is defined by the ``--snp-context-bases`` argument (specified individually for single base SNVs and indels; defaults 10 and 30 respectively).

Next the neural network output corresponding to the extracted sequence is extracted.
The fuzzy reference anchoring described above identifies the range of the neural network output containing the sequence of interest.

The sequence scoring function performs the forward algorithm over the neural network output to produce a score for the reference and proposed alternative sequence.
The difference between these two scores is the assigned score for the proposed variant.
Lower score are evidence for the alternative sequence and positive scores are evidence for the reference sequence.

These raw scores are softmax values over potential states, to match characteristics of a probability distribution.
In practice, these scores do not match emperical probabilities for a SNP given a truth dataset.
Thus a calibration step is applied to convert these scores to estimated emperical probabilities.
This allows more accurate estimation during aggregation across muliple reads at a genomic location.

Finally, calls across reads at each reference location are aggregated in order make a sampe reference call.
These results will be output into a VCF format file.

Currently ``diploid`` (default) and ``haploid`` SNP aggregation modes are available.
In ``haploid`` mode the probability of the reference and alternative alleles are simply the normalized product of the individual read probabilities.
In ``diploid`` mode the probability of each genotype (homozygous reference, heterozygous and homozygous alternative) are computed.
The probabilities for homozygous alleles are as in the ``haploid`` mode, while the heterozygous probability is given by the sum of the maximal probabilities taken over the sampling probability given a true diploid heterozygous allele.
These probabilities are then normalized given by Bayes' theorem.

---------------------
Modified Base Calling
---------------------

Modified base calling is performed largely in the same manner as SNP calling in terms of regions sequence associated nueral network output extraction.
The main difference is that instead of proposing alternative bases in the sequence, a modification is proposed.
This means that in order to identify a particular modification the model being used must be aware of this modification.
Training models for particular modifications of interest is described in the taiyaki documentation.

Use the ``--mod-motif`` argument in order to restrict tested locations to certain relevant motifs (e.g. ``--mod-motif Z CG 0`` to test only in CpG locations).
