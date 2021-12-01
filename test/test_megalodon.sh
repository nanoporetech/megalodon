############################
# Test files/dirs/settings #
############################

GUPPY_BIN_PATH="./ont-guppy-cpu/bin/"
GUPPY_FAST_CONFIG="dna_r9.4.1_450bps_fast.cfg"
GUPPY_MOD_CONFIG="dna_r9.4.1_450bps_modbases_5mc_hac.cfg"
MOD_CALIBRATION_FN="megalodon/megalodon/model_data/dna_r9.4.1_450bps_modbases_5mc_hac.cfg/megalodon_mod_calibration.npz"
TAIYAKI_CPT="taiyaki_model.checkpoint"
REMORA_MODEL="remora_model.onnx"
REMORA_SPEC="dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0"

FASTA_REF="reference.fa"
MINIMAP_INDEX="reference.fa.mmi"

VARS="variants.vcf.gz"

CTRL_READS="amplified_reads"
NAT_READS="native_reads"

NPROC=4


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
    --outputs basecalls \
    `# guppy options` \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_FAST_CONFIG} \
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
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_FAST_CONFIG} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# mapping settings (cram requires FASTA reference)` \
    --sort-mappings \
    --mappings-format cram \
    --cram-reference ${FASTA_REF} \
    `# per-read reference/signal mapping settings` \
    --ref-signal-mapping-offset 1 \
    --ref-length-range 500 3000 \
    --ref-percent-identity-threshold 90 \
    --ref-percent-coverage-threshold 90

# test mods output
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.mods_output \
    --overwrite \
    `# output all the mod things` \
    --outputs per_read_mods mods mod_mappings \
    `# guppy options` \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# modified base settings` \
    --mod-motif m CCWGG 1 \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mods-text \
    --write-mod-log-probs

# test remora mods output
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.remora_mods_output \
    --overwrite \
    `# output all mod the things` \
    --outputs per_read_mods mods mod_mappings mod_basecalls \
    `# guppy options` \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_FAST_CONFIG} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# modified base settings` \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mods-text \
    --write-mod-log-probs \
    `# remora settings` \
    --remora-model ${REMORA_MODEL}

# test example command
megalodon \
    ${CTRL_READS} \
    --guppy-config ${GUPPY_FAST_CONFIG} \
    --remora-modified-bases ${REMORA_SPEC} \
    --outputs basecalls mappings mod_mappings mods \
    --reference ${MINIMAP_INDEX} \
    --processes ${NPROC} \
    --overwrite

# test remora auto-load mods output
megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.remora_mods_output \
    --overwrite \
    `# output all mod the things` \
    --outputs per_read_mods mods mod_mappings mod_basecalls \
    `# guppy options` \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_FAST_CONFIG} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# modified base settings` \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mods-text \
    --write-mod-log-probs \
    `# remora settings` \
    --remora-modified-bases ${REMORA_SPEC}

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
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# mapping settings (cram requires FASTA reference)` \
    --sort-mappings \
    --mappings-format cram \
    --cram-reference ${FASTA_REF} \
    `# modified base settings` \
    --mod-motif m CCWGG 1 \
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
    --ref-percent-coverage-threshold 90

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
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    --guppy-timeout 300 \
    `# number of megalodon workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# mapping settings (cram requires FASTA reference)` \
    --sort-mappings \
    --mappings-format cram \
    --cram-reference ${FASTA_REF} \
    `# modified base settings` \
    --mod-motif m CCWGG 1 \
    --mod-output-formats bedmethyl modvcf wiggle \
    --write-mods-text \
    --write-mod-log-probs \
    `# sequence variant settings` \
    --variant-filename ${VARS} \
    --haploid \
    --write-variants-text \
    --write-vcf-log-probs \
    `# per-read reference/signal mapping settings` \
    --ref-length-range 500 30000 \
    --ref-percent-identity-threshold 90 \
    --ref-percent-coverage-threshold 90 \
    --ref-mods-all-motifs m 5mC CCWGG 1


##########################
# megalodon_extras tests #
##########################

megalodon_extras \
    modified_bases describe_alphabet \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_MOD_CONFIG}

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
    modified_bases split_by_motif ${FASTA_REF} \
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
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server
megalodon_extras \
    calibrate variants \
    --processes ${NPROC} \
    --out-pdf megalodon_variant_calibration.pdf \
    --overwrite

megalodon_extras \
    per_read_text modified_bases \
    ${CTRL_READS}.mega_res \
    --out-filename \
    ${CTRL_READS}.mega_res/per_read_modified_base_calls.rerun.txt
megalodon_extras \
    per_read_text variants \
    ${CTRL_READS}.mega_res \
    --out-filename \
    ${CTRL_READS}.mega_res/per_read_variant_calls.rerun.txt

sort -k1V -k2n ${NAT_READS}.mega_res/modified_bases.5mC.bed > \
     ${NAT_READS}.mega_res/modified_bases.5mC.sorted.bed
megalodon_extras \
    modified_bases per_site_thresholds \
    ${NAT_READS}.mega_res/ \
    ${NAT_READS}.mega_res/modified_bases.5mC.sorted.bed \
    --mod-bases m \
    --ground-truth-cov-min 3 --nanopore-cov-min 5

# test model calibration from mapped signal files
megalodon_extras \
    calibrate generate_mod_stats_from_msf \
    ${CTRL_READS}.mega_res/signal_mappings.hdf5 \
    --motif CCWGG 1 \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --out-filename ctrl_mod_stats.npz --modified-bases-set m \
    --processes 2
megalodon_extras \
    calibrate generate_mod_stats_from_msf \
    ${NAT_READS}.mega_res/signal_mappings.hdf5 \
    --motif CCWGG 1 \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --out-filename nat_mod_stats.npz --modified-bases-set m \
    --processes 2
megalodon_extras \
    calibrate merge_modified_bases_stats \
    ctrl_mod_stats.npz nat_mod_stats.npz \
    --out-filename mod_stats.all.npz
megalodon_extras \
    calibrate modified_bases \
    --ground-truth-llrs mod_stats.all.npz \
    --out-filename mod_calib.all.npz \
    --out-pdf mod_calib.all.pdf \
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
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server \
    --guppy-config ${GUPPY_MOD_CONFIG} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# modified base settings` \
    --mod-motif m CCWGG 1 \
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


#######################
# test other backends #
#######################

megalodon \
    `# input reads` \
    ${CTRL_READS} \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.taiyaki \
    --overwrite \
    `# output all the things` \
    --outputs basecalls mappings signal_mappings \
    `# taiyaki options` \
    --do-not-use-guppy-server \
    --taiyaki-model-filename ${TAIYAKI_CPT} \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX}

${GUPPY_BIN_PATH}/guppy_basecaller \
    -i ${CTRL_READS} \
    -s ${CTRL_READS}.post_out \
    -c ${GUPPY_MOD_CONFIG} \
    --fast5_out --post_out
megalodon \
    `# input reads` \
    ${CTRL_READS}.post_out \
    `# output location + overwrite` \
    --output-directory ${CTRL_READS}.post_out_mega_res \
    --overwrite \
    `# output all the things` \
    --outputs basecalls mappings \
    per_read_mods mods mod_mappings \
    `# FAST5 post out options` \
    --do-not-use-guppy-server \
    `# number of megalodon read processing workers` \
    --processes ${NPROC} \
    `# minimap2 index reference (recommended for memory efficiency)` \
    --reference ${MINIMAP_INDEX} \
    `# modified base settings` \
    --mod-motif m CCWGG 1 \
    --write-mods-text \
    --mod-calibration-filename \
    ${MOD_CALIBRATION_FN}


############
# test API #
############

python \
    megalodon/test/test_api.py \
    --guppy-server-path ${GUPPY_BIN_PATH}/guppy_basecall_server
