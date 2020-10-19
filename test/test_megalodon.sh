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
GUPPY_TIMEOUT=240


######################
# Main command tests #
######################

# test basecalls output
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.basecalls_output \
    --overwrite \
    `# output just basecalls` \
    --outputs basecalls mod_basecalls \
    `# guppy options` \
    --guppy-server-path ${GUPPY_PATH} \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    `# number of megalodon workers` \
    --processes ${NPROC}

# test mappings output
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.mappings_output \
    --overwrite \
    `# output all the things` \
    --outputs mappings per_read_refs signal_mappings \
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
    `# per-read reference/signal mapping settings` \
    --ref-length-range 500 3000 \
    --ref-percent-identity-threshold 90 \
    --ref-percent-coverage-threshold 90 \
    --ref-mods-all-motifs m 5mC CCWGG 1 \
    --ref-mods-all-motifs a 6mA GATC 1

# test mods output
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.mods_output \
    --overwrite \
    `# output all the things` \
    --outputs per_read_mods mods mod_mappings \
    `# guppy options` \
    --guppy-server-path ${GUPPY_PATH} \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    --guppy-timeout ${GUPPY_TIMEOUT} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# modified base settings` \
    --mod-motif Z CCWGG 1 \
    --mod-motif Y GATC 1 \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mods-text \
    --write-mod-log-probs

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

megalodon_extras \
    calibrate generate_variant_stats \
    ${CTRL_READS} \
    --reference ${MINIMAP_INDEX} \
    --num-reads 10 \
    --processes ${NPROC} \
    --guppy-server-path ${GUPPY_PATH}
megalodon_extras \
    calibrate variants \
    --processes ${NPROC} \
    --out-pdf megalodon_variant_calibration.pdf \
    --overwrite

# TODO add tests for more megalodon_extras commands


##############################
# test megalodon as pipeline #
##############################

# test running megalodon in separate steps
#   - useful to minimize time on GPU compute resources
megalodon \
    `# input reads (limit reads for tiny test)` \
    ${CTRL_READS} \
    --num-reads 5 \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.pipeline \
    --overwrite \
    `# output only per_read databases` \
    --outputs per_read_mods per_read_variants \
    `# guppy options` \
    --guppy-server-path ${GUPPY_PATH} \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    --guppy-timeout ${GUPPY_TIMEOUT} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# modified base settings` \
    --mod-motif Z CCWGG 1 \
    --mod-motif Y GATC 1 \
    `# sequence variant settings` \
    --variant-filename ${VARS} \
    `# skip database index to run as a pipeline` \
    --skip-database-index
# Create modified base database index
#   - Can be performed on CPU-only compute resources
megalodon_extras \
    modified_bases index_database \
    --megalodon-directory ${CTRL_READS}.pipeline
# Create sequence variants database index (not currently implemented)
#megalodon_extras \
#    variants index_database \
#    --megalodon-directory ${CTRL_READS}.pipeline
# aggregate mods and variants
#   - Can be performed on CPU-only compute resources
megalodon_extras \
    aggregate run \
    `# specify output options` \
    --megalodon-directory ${CTRL_READS}.pipeline \
    --output-suffix pipeline \
    --outputs mods variants \
    `# compute resources` \
    --processes ${NPROC} \
    `# modified base ouput options` \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mod-log-probs \
    `# sequence variant ouput options` \
    --haploid
