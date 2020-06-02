#!/bin/bash

TEST_DATA_DIR=${PWD}/data
echo "Looking for test data in ${TEST_DATA_DIR}"

#ORIG_DATA_FILE="${TEST_DATA_DIR}/ica_donor_5_channel_1.h5ad"
ORIG_DATA_FILE="${TEST_DATA_DIR}/ica_cord_blood_full_data.h5ad"
ORIG_DATA_BASENAME="${TEST_DATA_DIR}/ica_cord_blood_full_data"


CMD="bash ../../module/run_module.sh \
    --python.executable=python3 \
    --src.path=../../module/ \
    --data.file=${ORIG_DATA_FILE} \
    --output.basename=${ORIG_DATA_BASENAME} \
    --annotate=1 \
    --cells.min.counts=1000 \
    --cells.max.counts=20000 \
    --cells.min.genes=500 \
    --genes.min.cells=10 \
    --normalize=1 \
    --n.high.variance.genes=3000 \
    --compute.umap=1 \
    --compute.tsne=1 \
"

eval ${CMD}
