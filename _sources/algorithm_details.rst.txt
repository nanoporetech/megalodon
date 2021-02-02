***************************
Megalodon Algorithm Details
***************************

This page describes the details of how megalodon processes the raw nanopore signal to produce highly-accurate modified base and sequence variant calls.

------------
Base Calling
------------

Basecalling is performed exactly as in guppy.
Raw nanopore signal is normalized, chunked, processed with a recurrent neural network and decoded using Viterbi decoding.
Currently megalodon is only compatible with flip-flop basecalling networks (excluding RLE and k-mer based networks)
See `guppy documentation on the community page (login required) <https://community.nanoporetech.com/protocols/Guppy-protocol>`_ for more details.

-------------------
Reference Anchoring
-------------------

Megalodon's functionality centers on the anchoring of the high-information neural network basecalling output to a reference sequence.
Given anchored neural network output, alternatives to the reference (either modified bases or canonical bases) are proposed and scored to produce the highest accuracy results.

The neural network output is anchored to the reference via standard read mapping of produced basecalls to the reference sequence (maintaining the link to the neural network outputs).
If no reference mapping is produced (using ``minimap2`` via the ``mappy`` python interface) that read is not processed further (basecalls will be output if requested).
This standard read mapping is processed to produce a matching of each basecall with a reference position.
Reference positions within an insertion or deletion are assigned to the previous mapped read position (left justified; this behavior may change in future versions).
This constitutes the reference anchoring used for modified base and sequence variant calling steps.

------------------------
Sequence Variant Calling
------------------------

Megalodon currently filters alleles over a certain maximum size (default ``50``) as performance on larger indels has not currently been validated.
Note also that variants are converted into an "atomic" form (containing minimal unique variant sequence for indels).
Thus atomic variants do not contain context sequence and are expanded to include regions of ambiguity (indel within a repetitive region).

At each valid variant a region of context sequence around the variant is extracted.
The context sequence allows the scoring algorithm to traverse slightly different paths through the local neural network output.
The width of this sequence of interest is defined by the ``--variant-context-bases`` argument (specified individually for single base and insertion/deletion variants; defaults ``10`` and ``30`` respectively).

Next the neural network output corresponding to the reference sequence of interest is extracted.
The fuzzy reference anchoring described above identifies the range of the neural network output containing the sequence of interest.

The sequence scoring function performs the forward-backward algorithm and Viterbi decoding over the neural network output to produce a score for the reference and proposed alternative sequence.
The difference between these two scores is the assigned score for the proposed variant.
Lower (negative) score are evidence for the alternative sequence and higher (positive) scores are evidence for the reference sequence.

These raw scores are softmax values over potential states, to match characteristics of a probability distribution.
In practice, these scores do not match empirical probabilities for a variant given a truth data set.
Thus a calibration step is applied to convert these scores to estimated empirical probabilities.
This enables more accurate aggregation across reads.

As of version 1.0.0, megalodon now performs a second round of variant detection taking nearby variants into account.
Variants from the first round (considering each variant in isolation) are filtered to by a minimal probability of evidence for variant allele (default ``0.05``; set with ``--context-min-alt-prob`` argument).
In the second pass, variants within a set region are considered when estimating the probability of a particular variant (up to a set maximum number of context variants in order to reduce compute).
Scores for each potential context are combined statistically (using logsumexp) and these are the final scores reported for each variant.
This process reduces the number of false positives where a true variant is adjacent to another proposed variant.

Finally, calls across reads at each reference location are aggregated in order make a sample-level call.
These results will be output into a VCF format file.

Currently ``diploid`` (default) and ``haploid`` variant aggregation modes are available.
In ``haploid`` mode the probability of the reference and alternative alleles are simply the normalized (via Bayes' theorem) product of the individual read probabilities.
In ``diploid`` mode the probability of each genotype (homozygous reference, heterozygous and homozygous alternative) are computed.
The probabilities for homozygous alleles are as in the ``haploid`` mode, while the heterozygous probability is given by the weighted sum of the maximal probabilities taken over the sampling distribution (binomial with ``p=0.5``) given a true diploid heterozygous allele.
These probabilities are then normalized given by Bayes' theorem.

---------------------
Modified Base Calling
---------------------

Modified base calling is performed largely in the same manner as variant calling above in terms of sequence and associated neural network output extraction.
The main difference is that instead of proposing alternative bases in the sequence, a modification is proposed.
This means that in order to identify a particular modification the model must be aware of this modification.
Training models for particular modifications of interest is described in `Megalodon documentation here <https://nanoporetech.github.io/megalodon/model_training.html>`_.

Use the ``--mod-motif`` argument in order to restrict tested locations to certain relevant motifs (e.g. ``--mod-motif Z CG 0`` to test only in CpG locations).

Modified bases can also be output anchored to the basecalls as in Guppy, but these calls are generally not as accurate than the reference anchored calls.
These ``mod_basecalls`` are output in the BAM ``Mm`` and ``Ml`` tags as specified by hts-specs.
