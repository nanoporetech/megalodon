***********************************
``megalodon_extras phase_variants``
***********************************

The ``megalodon_extras phase_variants`` command group contains commands to assist in the Megalodon pipeline to produce the highest quality phased variants calls.

---------------------------------------------------
``megalodon_extras phase_variants whatshap_filter``
---------------------------------------------------

`WhatsHap <https://whatshap.readthedocs.io/en/latest/>`_ (as of version ``0.18``) cannot process some complex variants.
Providing such variants causes WhatsHap to exit with an error.
This command is provided to remove these complex variants and allow processing to proceed without error.
Note that these variants are still considered outside of the WhatsHap phasing step of the Megalodon phasing pipeline.

-----------------------------------------------------------
``megalodon_extras phase_variants extract_haplotype_reads``
-----------------------------------------------------------

From alignment files produced by ``whatshap haplotag``, extract read ids for reads assigned to one of the two haplotypes.
One file will be produced for each haplotype value in the alignment file (two for standard diploid processing).

----------------------------------------------------------
``megalodon_extras phase_variants merge_haploid_variants``
----------------------------------------------------------

Merge haploid calls from original Megalodon variants and separate haplotype sets of calls.

This command should only be used as recommended in the Megalodon variant phasing pipeline.
Use of this command outside of this context is not recommended as several processing steps depend upon the variant processing steps.

This command iterates over the three sets of sorted variants.
If a variant is not found in the haplotype variant files, the call from the original Megalodon run is taken.
If a variant is called as homozygous in the original Megalodon calling a heterozygous call cannot be output.
If a variant was originally calls as heterozygous and the variant is called in both haplotype calls, then the output call is determined from the two haplotype calls.
