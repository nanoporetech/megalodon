***************
Variant Phasing
***************

This page walks through the steps to use megalodon in conjunction with `whatshap <https://whatshap.readthedocs.io/en/latest/>`_ to produce the highest quality phased variant calls.

This workflow requires a working installation of whatshap. See instructions at the link above.

--------
Workflow
--------

::

   # run megalodon to produce whatshap_mappings
   megalodon \
       fast5s --outputs whatshap_mappings --processes 4 \
       --reference reference.fasta --variant-filename variants.vcf.gz --overwrite

   # run whatshap with produced mappings and variants
   samtools index megalodon_results/whatshap_mappings.sorted.bam
   whatshap \
       phase --indels --distrust-genotypes \
       -o megalodon_results/variants.phased.vcf \
       megalodon_results/variants.sorted.vcf.gz \
       megalodon_results/whatshap_mappings.sorted.bam

   # color reads against the phased variants
   bgzip megalodon_results/variants.phased.vcf
   tabix megalodon_results/variants.phased.vcf.gz
   whatshap \
       haplotag megalodon_results/variants.phased.vcf.gz \
       megalodon_results/whatshap_mappings.sorted.bam \
       -o megalodon_results/whatshap_mappings.haplotagged.bam

   # extract relevant reads and call haploid variants (more accurate)
   python \
       megalodon/scripts/extract_haplotype_read_ids.py \
       megalodon_results/whatshap_mappings.haplotagged.bam \
       megalodon_results/whatshap_mappings.haploid_reads
   python \
       megalodon/scripts/run_aggregation.py \
       --outputs snps --haploid --output-suffix haplotype_1 \
       --read-ids-filename megalodon_results/whatshap_mappings.haploid_reads.haplotype_1_read_ids.txt \
       --reference reference.fasta
   python \
       megalodon/scripts/run_aggregation.py \
       --outputs snps --haploid --output-suffix haplotype_2 \
       --read-ids-filename megalodon_results/whatshap_mappings.haploid_reads.haplotype_2_read_ids.txt \
       --reference reference.fasta

   # merge haploid variants to produce diploid variants
   python \
       megalodon/scripts/merge_haploid_variants.py \
       megalodon_results/variants.sorted.vcf.gz \
       megalodon_results/variants.haplotype_1.sorted.vcf.gz \
       megalodon_results/variants.haplotype_2.sorted.vcf.gz \
       --out-vcf megalodon_results/variants.haploid_merged.vcf
