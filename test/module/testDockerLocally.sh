#!/bin/bash

TEST_DATA_DIR=${PWD}/data
echo "Looking for test data in ${TEST_DATA_DIR}"

docker build -t scanpy_utilities_module ../..

ORIG_DATA_FILE="${TEST_DATA_DIR}/ica_donor_5_channel_1.h5ad"
ORIG_DATA_BASENAME="${TEST_DATA_DIR}/ica_donor_5_channel_1"
eval "rm -f ${TEST_DATA_DIR}/ica_donor_5_channel_1_*"

CMD="bash build/run_module.sh \
    --src.path=build/ \
    --data.file=${ORIG_DATA_FILE} \
    --output.basename=${ORIG_DATA_BASENAME} \
    --annotate=1 \
    --cells.min.counts=1000 \
    --cells.max.counts=20000 \
    --cells.min.genes=500 \
    --genes.min.cells=10 \
    --cell.type.marker.file=${TEST_DATA_DIR}/markers.txt \
    --normalize=1 \
    --n.high.variance.genes=3000 \
    --compute.umap=1 \
    --compute.tsne=1 \
"

docker run -v ${TEST_DATA_DIR}:${TEST_DATA_DIR} -t scanpy_utilities_module ${CMD}
