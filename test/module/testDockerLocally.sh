#!/bin/bash

TEST_DATA_DIR=${PWD}/data
echo "Looking for test data in ${TEST_DATA_DIR}"

#docker build -t scanpy_utilities_module ../..

#ORIG_DATA_FILE="${TEST_DATA_DIR}/ica_donor_5_channel_1.h5ad"
#ORIG_DATA_BASENAME="${TEST_DATA_DIR}/ica_donor_5_channel_1"
ORIG_DATA_FILE="${TEST_DATA_DIR}/ica_cord_blood_full_data"
ORIG_DATA_BASENAME="${TEST_DATA_DIR}/output_full_data"


eval "rm -f ${TEST_DATA_DIR}/ica_donor_5_channel_1_*"

BASE_CMD="bash build/run_module.sh \
    --src.path=build/ \
    --output.basename=${ORIG_DATA_BASENAME} \
"

TEST_1="${BASE_CMD} \
    --data.file=${ORIG_DATA_FILE} \
    --annotate=1 \
    --cells.min.counts=1000 \
    --cells.max.counts=20000 \
    --cells.min.genes=500 \
    --genes.min.cells=10 \
    --cell.type.marker.file=${TEST_DATA_DIR}/markers.txt \
    --gene.annotation.database=org.Hs.eg.db \
    --normalize=1 \
    --n.high.variance.genes=3000 \
    --compute.umap=1 \
    --compute.tsne=1 \
"

#TEST_2="${BASE_CMD} \
#    --data.file=${ORIG_DATA_BASENAME}_gene_filter.h5ad \
#    --cell.type.marker.file=${TEST_DATA_DIR}/ensembl_markers.txt
#"

#TEST_3="${BASE_CMD} \
#    --data.file=${TEST_DATA_DIR}/ica_cord_blood_donor_5_filtered.h5ad \
#    --cell.type.marker.file=${TEST_DATA_DIR}/ica_cord_blood_markers.txt \
#    --gene.annotation.database=org.Hs.eg.db
#"


declare -a tests=(
    "${TEST_1}"
#    "${TEST_3}"
)

for tst in "${tests[@]}"
do
    echo "========"
    echo "NEW TEST"
    echo "========"
    echo docker run -v ${TEST_DATA_DIR}:${TEST_DATA_DIR} -t scanpy_utilities_module ${tst}
    docker run -v ${TEST_DATA_DIR}:${TEST_DATA_DIR} -t scanpy_utilities_module ${tst}
done
