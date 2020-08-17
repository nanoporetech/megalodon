*****************************
``megalodon_extras variants``
*****************************

The ``megalodon_extras variants`` command group contains various commands related to sequence variant processing within Megalodon.

-------------------------------------
``megalodon_extras variants atomize``
-------------------------------------

This command processes each variant in a VCF file to convert the variants to their atomic form.
When the variant file produced by this command is used in the ``megalodon`` command with the ``--variants-are-atomized`` flag, processing speed can be drastically increased.
This is especially true for high-density variant files.

For each variant processed in Megalodon the "atomic" variant is required for the highest quality results.
For single nucleotide variants (SNVs) the atomic form simply reduces any multi-nucleotide SNP into the single nucleotide swaps.
For insertions and deletions (indels) the atomic form removes context bases and expands the indel to include all unambiguous positions.
For example consider a bit of reference sequence ``AGCGCA``.
An insertion after the first base of a ``GC`` could be validly annotated in a VCF file as reference allele ``A`` and alternative allele ``AGC``.
But the atomic form of this variant would be reference allele ``GCGC`` and alternative allele ``GCGCGC``.
In this way atomic variants capture all reference bases impacted by a variant and no more.
Processing variants in this way ensures that the correct bit of sequence and raw signal is considered when determining the correct allele within Megalodon.

The processing of each variant to its atomic form must be completed each time a variant overlaps a read since VCF files can often be too large to store in RAM.
Thus the computational cost for high coverage or dense variant files can be quite high for this step of the processing.
This command allows this atomizing step to be performed once before starting a ``megalodon`` run.
When passing a variant file processed with this command, the ``--variants-are-atomized`` flag should be specified.
Note that this flag should not be used with VCF files originating from other sources.
Specifically, this command adds a non-standard ``CB`` (context base) VCF flag in order to indicate that a context base added to a variant is not part of the atomic form of that variant.

-------------------------------------
``megalodon_extras variants resolve``
-------------------------------------

The ``megalodon_extras variants resolve`` command resolves conflicting variants and removes variants determined to be the reference allele.

Megalodon processes sequence variants to make a call for each variant individually (taking nearby variants into account).
Nearby variants are often overlapping and thus incompatible.
The primary use for this command is to resolve situations where multiple overlapping variants are called.
At each site with overlapping variants, the one with the highest probability is selected for the output file.

In addition, this command has options to filter variants which were called as the reference allele (``--max-likelihood-ratio``), filter for coverage (``--min-depth``), and revert atomic variant notation (``--trim-variants``).

There are also a number of options to inspect potential systematic bias in identified sequence variants.
Providing a VCF file called only from reverse strand mapping reads via the ``--reverse-strand-variants`` argument activates the output.
The main VCF provided is then assumed to be derived from forward strand mapping reads only.
In this mode, variants are output when they are identified only on one strand and not the other to allow analysis of potential bias in basecalling models.
This feature is experimental and does not have defined pipelines for downstream use.

-------------------------------------------------
``megalodon_extras variants heterozygous_factor``
-------------------------------------------------

Determine the result of the heterozygous factor on identifying the correct balance of homozygous to heterozygous variants.

This command can assist in setting an optimal value for the ``--heterozygous-factors`` argument.
The default value is intended to minimize the number of false heterozygous calls.
The recommended phased variant pipeline for Megalodon processes variants such that if a variant is not initially called heterozygous, a heterozygous call cannot be made, but heterozygous calls can be converted to homozygous calls.
Since this phased variant pipeline is recommended in order to obtain the highest quality variants, false homozygous calls are minimized.

If the aim is to achieve a balance of homozygous and heterozygous calls, this command can be used to evaluate this balance for a particular ``--heterozygous-factors`` setting.
This command will output the number of called homozygous and heterozygous calls compared to their ground truth from a provided set of variants.
As a guide, previous Megalodon version had the default value of ``--heterozygous-factors 2.1 1.6`` which achieves a better balance than the current default of ``--heterozygous-factors 1.0 1.0`` which minimizes false homozygous calls.

--------------------------------------------
``megalodon_extras variants index_database``
--------------------------------------------

This command is not currently implemented, but will be in a future release.
