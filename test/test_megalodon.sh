############################
# Test files/dirs/settings #
############################

GUPPY_PATH="./ont-guppy-cpu/bin/guppy_basecall_server"
GUPPY_FAST_CONFIG="dna_r9.4.1_450bps_fast.cfg"
GUPPY_MOD_CONFIG="dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg"

FASTA_REF="reference.fa"
MINIMAP_INDEX="reference.fa.mmi"

VARS="variants.vcf.gz"

CTRL_READS="amplified_reads"
NAT_READS="native_reads"

NPROC=8
GUPPY_TIMEOUT=60


######################
# Main command tests #
######################

# test simple basecalling
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.basecalls_only.mega_res \
    --overwrite \
    `# output just basecalls` \
    --outputs basecalls \
    `# guppy options` \
    --guppy-server-path ${GUPPY_PATH} \
    --guppy-config ${GUPPY_FAST_CONFIG} \
    `# number of megalodon workers` \
    --processes ${NPROC}

# test outputting everything
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.mega_res \
    --overwrite \
    `# output all the things` \
    --outputs basecalls mod_basecalls mappings \
    per_read_mods mods mod_mappings \
    per_read_variants variants variant_mappings \
    per_read_refs signal_mappings \
    `# guppy options` \
    --guppy-server-path ${GUPPY_PATH} \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    --guppy-timeout ${GUPPY_TIMEOUT} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# mapping settings (cram requires FASTA reference)` \
    --sort-mappings \
    --mappings-format cram \
    --cram-reference ${FASTA_REF} \
    `# modified base settings` \
    --mod-motif Z CCWGG 1 \
    --mod-motif Y GATC 1 \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mods-text \
    --write-mod-log-probs \
    `# sequence variant settings` \
    --variant-filename ${VARS} \
    --haploid \
    --write-variants-text \
    --write-vcf-log-probs \
    `# per-read reference/signal mapping settings` \
    --ref-length-range 500 3000 \
    --ref-percent-identity-threshold 90 \
    --ref-percent-coverage-threshold 90 \
    --ref-mods-all-motifs m 5mC CCWGG 1 \
    --ref-mods-all-motifs a 6mA GATC 1

# process native reads for downstream results
megalodon \
    `# input reads` \
    ${NAT_READS} \
    `# output location + overwrite` \
    --output-directory ${NAT_READS}.mega_res \
    --overwrite \
    `# output all the things` \
    --outputs basecalls mod_basecalls mappings \
    per_read_mods mods mod_mappings \
    per_read_variants variants variant_mappings \
    per_read_refs signal_mappings \
    `# guppy options` \
    --guppy-server-path ${GUPPY_PATH} \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    --guppy-timeout ${GUPPY_TIMEOUT} \
    `# number of megalodon workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# mapping settings (cram requires FASTA reference)` \
    --sort-mappings \
    --mappings-format cram \
    --cram-reference ${FASTA_REF} \
    `# modified base settings` \
    --mod-motif Z CCWGG 1 \
    --mod-motif Y GATC 1 \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mods-text \
    --write-mod-log-probs \
    `# sequence variant settings` \
    --variant-filename ${VARS} \
    --haploid \
    --write-variants-text \
    --write-vcf-log-probs \
    `# per-read reference/signal mapping settings` \
    --ref-length-range 500 3000 \
    --ref-percent-identity-threshold 90 \
    --ref-percent-coverage-threshold 90 \
    --ref-mods-all-motifs m 5mC CCWGG 1 \
    --ref-mods-all-motifs a 6mA GATC 1


##########################
# megalodon_extras tests #
##########################

# add extras commands to test
megalodon_extras \
    aggregate run --outputs variants mods \
    --megalodon-directory ${CTRL_READS}.mega_res \
    --processes ${NPROC}

megalodon_extras \
    merge modified_bases ${CTRL_READS}.mega_res/ ${NAT_READS}.mega_res \
    --output-megalodon-results-dir merge_mods \
    --overwrite
megalodon_extras \
    aggregate run --outputs mods \
    --megalodon-directory merge_mods \
    --processes ${NPROC}

rm -r megalodon_results.split_by_motif.CCAGG_1 \
   megalodon_results.split_by_motif.CCTGG_1
megalodon_extras \
    modified_bases split_by_motif e_coli_reference.fa \
    --motif CCAGG 1 --motif CCTGG 1 \
    --megalodon-directory ${NAT_READS}.mega_res/
megalodon_extras \
    aggregate run \
    --megalodon-directory megalodon_results.split_by_motif.CCAGG_1 \
    --outputs mods \
    --processes ${NPROC}
megalodon_extras \
    aggregate run \
    --megalodon-directory megalodon_results.split_by_motif.CCTGG_1 \
    --outputs mods \
    --processes ${NPROC}

# TODO add tests for more megalodon_extras commands
