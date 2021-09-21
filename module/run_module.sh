#!/bin/bash

# version information
echo "Version v0.1"

# default values for parameters
PY_EXEC="python"
SRC_PATH=.
COMPUTE_UMAP=0
COMPUTE_TSNE=0

# helper function
function exitOnError(){
    if [ $1 -ne  0 ]; then
        echo "ERROR CODE $1"
        echo $2
        exit $1
    fi
}

# parse command line parameters
for i in "$@"
do
    case $i in
        --python.executable=*)
        PY_EXEC="${i#*=}"
        echo "Using ${PY_EXEC} for python executable"
        ;;

        --src.path=*)
        SRC_PATH="${i#*=}"
        echo "Looking for code in ${SRC_PATH}"
        ;;

        --data.file=*)
        DATA_FILE="${i#*=}"
        echo "data file: ${DATA_FILE}"
        ;;

        --mito.file=*)
        MITO_FILE="${i#*=}"
        echo "mitocondrial genes file: ${MITO_FILE}"
        ;;

        --output.basename=*)
        OUTPUT_BASENAME="${i#*=}"
        echo "output basename: ${OUTPUT_BASENAME}"
        ;;

        --genome=*)
        GENOME="${i#*=}"
        echo "genome: ${GENOME}"
        ;;

        --annotate=*)
        ANNOTATE="${i#*=}"
        echo "annotate: ${ANNOTATE}"
        ;;

        --cells.min.counts=*)
        CELLS_MIN_COUNTS="${i#*=}"
        echo "cells min counts: ${CELLS_MIN_COUNTS}"
        ;;

        --cells.max.counts=*)
        CELLS_MAX_COUNTS="${i#*=}"
        echo "cells max counts: ${CELLS_MAX_COUNTS}"
        ;;

        --cells.min.genes=*)
        CELLS_MIN_GENES="${i#*=}"
        echo "cells min genes: ${CELLS_MIN_GENES}"
        ;;

        --cells.max.genes=*)
        CELLS_MAX_GENES="${i#*=}"
        echo "cells max genes: ${CELLS_MAX_GENES}"
        ;;

        --cells.max.mt.pct=*)
        CELLS_MAX_MT_PCT="${i#*=}"
        echo "cells max mitocondrial percent: ${CELLS_MAX_MT_PCT}"
        ;;

        --genes.min.counts=*)
        GENES_MIN_COUNTS="${i#*=}"
        echo "genes min counts: ${GENES_MIN_COUNTS}"
        ;;

        --genes.max.counts=*)
        GENES_MAX_COUNTS="${i#*=}"
        echo "genes max counts: ${GENES_MAX_COUNTS}"
        ;;

        --genes.min.cells=*)
        GENES_MIN_CELLS="${i#*=}"
        echo "genes min cells: ${GENES_MIN_CELLS}"
        ;;

        --genes.max.cells=*)
        GENES_MAX_CELLS="${i#*=}"
        echo "genes max cells: ${GENES_MAX_CELLS}"
        ;;

        --cell.type.marker.file=*)
        CELL_TYPE_MARKER_FILE="${i#*=}"
        echo "cell type marker file: ${CELL_TYPE_MARKER_FILE}"
        ;;

        --gene.annotation.database=*)
        GENE_ANNO_DB="${i#*=}"
        echo "gene annotation database: ${GENE_ANNO_DB}"
        ;;
    
        --normalize=*)
        NORMALIZE="${i#*=}"
        echo "normalize: ${NORMALIZE}"
        ;;

        --batch.correct=*)
        BATCH_CORRECT="${i#*=}"
        echo "batch correct: ${BATCH_CORRECT}"
        ;;

        --batch.variable=*)
        BATCH_VAR="${i#*=}"
        echo "batch variable: ${BATCH_VAR}"
        ;;

        --n.high.variance.genes=*)
        HIGH_VAR_GENES="${i#*=}"
        echo "num high variance genes: ${HIGH_VAR_GENES}"
        ;;

        --compute.umap=*)
        COMPUTE_UMAP="${i#*=}"
        echo "compute umap: ${COMPUTE_UMAP}"
        ;;

        --compute.tsne=*)
        COMPUTE_TSNE="${i#*=}"
        echo "compute tsne: ${COMPUTE_TSNE}"
        ;;

        *)
            echo "ERROR: unrecognized option"
            echo "${i%=*}"
            exit 1
        ;;
    esac
done

if [ ! -z ${GENOME} ] && [ ${DATA_FILE: -3} == ".h5" ]; then
    echo "-- converting from H5 to H5AD format --"
    FULL_OUTPUT="${OUTPUT_BASENAME}.h5ad"
    eval "${PY_EXEC} ${SRC_PATH}/convert_to_h5ad.py ${DATA_FILE} ${GENOME} ${FULL_OUTPUT}"
    exitOnError $? "Error during conversion from h5 to h5ad format."
    DATA_FILE=${FULL_OUTPUT}
fi

if [ ! -z ${ANNOTATE} ] && [ "${ANNOTATE}" -eq 1 ]; then
    echo "-- adding default annotations --"
    FULL_OUTPUT=${OUTPUT_BASENAME}"_annotated.h5ad"
    eval "${PY_EXEC} ${SRC_PATH}/add_default_annotations.py ${DATA_FILE} ${FULL_OUTPUT}"
    exitOnError $? "Error during adding default annotations."
    DATA_FILE=${FULL_OUTPUT}
fi

if [ ! -z ${CELLS_MIN_COUNTS} ] || [ ! -z ${CELLS_MAX_COUNTS} ] || [ ! -z ${CELLS_MIN_GENES} ] || [ ! -z ${CELLS_MAX_GENES} || [ ! -z ${MITO_FILE}] ; then
    echo "-- filtering cells --"
    if [ -z ${CELLS_MIN_COUNTS} ]; then CELLS_MIN_COUNTS="0"; fi
    if [ -z ${CELLS_MAX_COUNTS} ]; then CELLS_MAX_COUNTS="0"; fi
    if [ -z ${CELLS_MIN_GENES} ]; then CELLS_MIN_GENES="0"; fi
    if [ -z ${CELLS_MAX_GENES} ]; then CELLS_MAX_GENES="0"; fi
    if [ -z ${MITO_FILE} ]; then MITO_FILE="SKIP"; fi
    if [ -z ${CELLS_MAX_MT_PCT} ]; then CELLS_MAX_MT_PCT="100"; fi
    FULL_OUTPUT="${OUTPUT_BASENAME}_cell_filter.h5ad"
    eval "${PY_EXEC} ${SRC_PATH}/filter_cells.py ${DATA_FILE} ${FULL_OUTPUT} ${CELLS_MIN_COUNTS} ${CELLS_MAX_COUNTS} ${CELLS_MIN_GENES} ${CELLS_MAX_GENES} ${MITO_FILE} ${CELLS_MAX_MT_PCT}"
    exitOnError $? "Error during cell filtering."
    DATA_FILE=${FULL_OUTPUT}
fi

if [ ! -z ${GENES_MIN_COUNTS} ] || [ ! -z ${GENES_MAX_COUNTS} ] || [ ! -z ${GENES_MIN_CELLS} ] || [ ! -z ${GENES_MAX_CELLS} ] ; then
    echo "-- filtering genes --"
    if [ -z ${GENES_MIN_COUNTS} ]; then GENES_MIN_COUNTS="0"; fi
    if [ -z ${GENES_MAX_COUNTS} ]; then GENES_MAX_COUNTS="0"; fi
    if [ -z ${GENES_MIN_CELLS} ]; then GENES_MIN_CELLS="0"; fi
    if [ -z ${GENES_MAX_CELLS} ]; then GENES_MAX_CELLS="0"; fi
    FULL_OUTPUT="${OUTPUT_BASENAME}_gene_filter.h5ad"
    eval "${PY_EXEC} ${SRC_PATH}/filter_genes.py ${DATA_FILE} ${FULL_OUTPUT} ${GENES_MIN_COUNTS} ${GENES_MAX_COUNTS} ${GENES_MIN_CELLS} ${GENES_MAX_CELLS}"
    exitOnError $? "Error during gene filtering."
    DATA_FILE=${FULL_OUTPUT}
fi

if [ ! -z ${CELL_TYPE_MARKER_FILE} ]; then
    echo "-- labelling cell types --"
    FULL_OUTPUT="${OUTPUT_BASENAME}_cell_types.h5ad"
    echo "cell type indentification input: ${DATA_FILE}"
    echo "cell type indentification output: ${FULL_OUTPUT}"
    echo "cell type marker file: ${CELL_TYPE_MARKER_FILE}"
    cp ${DATA_FILE} ${FULL_OUTPUT}
    Rscript ${SRC_PATH}/identify_cell_types.R ${FULL_OUTPUT} ${CELL_TYPE_MARKER_FILE} ${GENE_ANNO_DB}
    exitOnError $? "Error during cell type identification"
    DATA_FILE=${FULL_OUTPUT}    
fi

if [ ! -z ${NORMALIZE} ] && [ "${NORMALIZE}" -eq 1 ]; then
    echo "-- normalizing --"
    FULL_OUTPUT="${OUTPUT_BASENAME}_normalized.h5ad"
    eval "${PY_EXEC} ${SRC_PATH}/generate_clusters_for_normalization.py ${DATA_FILE}"
    exitOnError $? "Error generating clusters for normalization."
    Rscript ${SRC_PATH}/compute_size_factors.R
    exitOnError $? "Error computing size factors for normalization."
    eval "${PY_EXEC} ${SRC_PATH}/normalize.py ${DATA_FILE} ${FULL_OUTPUT}"
    exitOnError $? "Error normalizing the data."
    rm -f temp_clustered_for_scran.h5ad
    rm -f temp_size_factors.csv
    DATA_FILE=${FULL_OUTPUT}
fi

if [ ! -z ${BATCH_CORRECT} ] && [ "${BATCH_CORRECT}" -eq 1 ]; then
    echo "-- batch correcting --"
    FULL_OUTPUT="${OUTPUT_BASENAME}_batch_corrected.h5ad"
    echo "batch correction input: ${DATA_FILE}"
    echo "batch correction output: ${FULL_OUTPUT}"
    echo "batch variable: ${BATCH_VAR}"
    cp ${DATA_FILE} ${FULL_OUTPUT}
    Rscript ${SRC_PATH}/batch_correct.R ${FULL_OUTPUT} ${BATCH_VAR}
    exitOnError $? "Error during batch correction."
    DATA_FILE=${FULL_OUTPUT}
fi

if [ ! -z ${HIGH_VAR_GENES} ]; then
    echo "-- selecting ${HIGH_VAR_GENES} high variance genes --"
    eval "${PY_EXEC} ${SRC_PATH}/high_variance_genes.py ${DATA_FILE} ${OUTPUT_BASENAME} ${HIGH_VAR_GENES}"
    exitOnError $? "Error filtering high variance genes."
    DATA_FILE="${OUTPUT_BASENAME}_high_variance_genes_annotated.h5ad"
fi

if [[ "${COMPUTE_UMAP}" -eq 1 || "${COMPUTE_TSNE}" -eq 1 ]]; then
    FULL_OUTPUT="${OUTPUT_BASENAME}_dim_reduce.h5ad"
    eval "${PY_EXEC} ${SRC_PATH}/dimension_reduction.py ${DATA_FILE} ${FULL_OUTPUT} ${COMPUTE_UMAP} ${COMPUTE_TSNE}"
    exitOnError $? "Error creating umap or tsne plots."
    DATA_FILE=${FULL_OUTPUT}
fi
