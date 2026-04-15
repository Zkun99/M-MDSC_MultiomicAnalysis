#default value
input_loom=inhouse_late_MDSC.loom
n_workers=1


# 数据库文件目录（根据你的实际情况修改）
db_dir="/data/0_gzk/sc_reference/transcription_factor_scenic/"

# 数据库文件名
tfs="allTFs_hg38_correct.txt"
feather="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
tbl="motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# 使用完整路径
tfs="${db_dir}${tfs}"
feather="${db_dir}${feather}"
tbl="${db_dir}${tbl}"

# PySCENIC路径（如果pyscenic已在PATH中，可以直接用"pyscenic"）
pyscenic="pyscenic"

# 输出目录（可选：创建输出目录）
output_dir="pyscenic_output"
mkdir -p "${output_dir}"

# ============================================
# 执行部分
# ============================================

echo "========================================"
echo "PySCENIC Analysis Start"
echo "========================================"
echo "Input loom: ${input_loom}"
echo "Workers: ${n_workers}"
echo "Database dir: ${db_dir}"
echo "Output dir: ${output_dir}"
echo "========================================"

# 检查文件是否存在
echo "Checking files..."
[ -f "${input_loom}" ] || { echo "Error: Input loom file not found: ${input_loom}"; exit 1; }
[ -f "${tfs}" ] || { echo "Error: TFs file not found: ${tfs}"; exit 1; }
[ -f "${feather}" ] || { echo "Error: Feather file not found: ${feather}"; exit 1; }
[ -f "${tbl}" ] || { echo "Error: Motif file not found: ${tbl}"; exit 1; }
echo "All required files found."

# 步骤1: GRN
echo ""
echo "Step 1: Running GRNBoost2..."
grn_output="${output_dir}/adjacencies.tsv"
${pyscenic} grn \
    --num_workers ${n_workers} \
    --output "${grn_output}" \
    --method grnboost2 \
    --client_or_address "local"\
    "${input_loom}" "${tfs}"
[ $? -eq 0 ] && [ -f "${grn_output}" ] || { echo "Error: GRNBoost2 failed"; exit 1; }
echo "GRNBoost2 completed: ${grn_output}"

# 步骤2: cisTarget
echo ""
echo "Step 2: Running cisTarget..."
ctx_output="${output_dir}/regulons.csv"
${pyscenic} ctx \
    "${grn_output}" "${feather}" \
    --annotations_fname "${tbl}" \
    --expression_mtx_fname "${input_loom}" \
    --mode "dask_multiprocessing" \
    --output "${ctx_output}" \
    --num_workers ${n_workers} \
    --mask_dropouts
[ $? -eq 0 ] && [ -f "${ctx_output}" ] || { echo "Error: cisTarget failed"; exit 1; }
echo "cisTarget completed: ${ctx_output}"

# 步骤3: AUCell
echo ""
echo "Step 3: Running AUCell..."
auc_output="${output_dir}/auc_mtx.csv"
${pyscenic} aucell \
    "${input_loom}" "${ctx_output}" \
    --output "${auc_output}" \
    --num_workers ${n_workers}
[ $? -eq 0 ] && [ -f "${auc_output}" ] || { echo "Error: AUCell failed"; exit 1; }
echo "AUCell completed: ${auc_output}"

echo ""
echo "========================================"
echo "PySCENIC Analysis Completed Successfully!"
echo "========================================"
echo "Output files in ${output_dir}:"
echo "  1. GRN network: ${grn_output}"
echo "  2. Regulons: ${ctx_output}"
echo "  3. AUC matrix: ${auc_output}"
echo "========================================"