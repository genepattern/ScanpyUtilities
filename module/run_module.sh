#!/bin/bash

echo "Version v0.1"
SRC_PATH=.
  
function exitOnError(){
    if [ $1 -ne  0 ]; then
        echo "ERROR CODE $1"
        echo $2
        exit $1
    fi
}


for i in "$@"
do
    case $i in
        --src.path=*)
        SRC_PATH="${i#*=}"
        echo "Looking for code in ${SRC_PATH}"
        ;;

        --data.file=*)
        DATA_FILE="${i#*=}"
        echo "data file: ${DATA_FILE}"
        ;;

        --output.basename=*)
        OUTPUT_BASENAME="${i#*=}"
        echo "output basename: ${OUTPUT_BASENAME}"
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

        --genome=*)
        GENOME="${i#*=}"
        echo "genome: ${GENOME}"
        ;;

        *)
            echo "ERROR: unrecognized option"
            echo "${i%=*}"
            exit 1
        ;;
    esac
done


if [ ! -z $GENOME ] && [ ${DATA_FILE: -3} == ".h5"   ]; then
    echo "-- converting from H5 to H5AD format --"
    FULL_OUTPUT=$OUTPUT_BASENAME".h5ad"
    python3 $SRC_PATH/convert_to_h5ad.py $DATA_FILE $GENOME $FULL_OUTPUT
    exitOnError $? "Error during conversion from h5 to h5ad format."
    DATA_FILE=$FULL_OUTPUT
fi


if [ ! -z $ANNOTATE ] && [ "$ANNOTATE" -eq 1 ]; then
    echo "-- adding default annotations --"
    FULL_OUTPUT=$OUTPUT_BASENAME"_annotated.h5ad"
    python3 $SRC_PATH/add_default_annotations.py $DATA_FILE $FULL_OUTPUT
    exitOnError $? "Error during adding default annotations."
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $CELLS_MIN_COUNTS ] || [ ! -z $CELLS_MAX_COUNTS ] || [ ! -z $CELLS_MIN_GENES ] || [ ! -z $CELLS_MAX_GENES ] ; then
    echo "-- filtering cells --"
    if [ -z $CELLS_MIN_COUNTS ]; then CELLS_MIN_COUNTS="0"; fi
    if [ -z $CELLS_MAX_COUNTS ]; then CELLS_MAX_COUNTS="0"; fi
    if [ -z $CELLS_MIN_GENES ]; then CELLS_MIN_GENES="0"; fi
    if [ -z $CELLS_MAX_GENES ]; then CELLS_MAX_GENES="0"; fi

    FULL_OUTPUT=$OUTPUT_BASENAME"_cell_filter.h5ad"
    python3 $SRC_PATH/filter_cells.py $DATA_FILE $FULL_OUTPUT $CELLS_MIN_COUNTS $CELLS_MAX_COUNTS $CELLS_MIN_GENES $CELLS_MAX_GENES
    exitOnError $? "Error during cell filtering."
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $GENES_MIN_COUNTS ] || [ ! -z $GENES_MAX_COUNTS ] || [ ! -z $GENES_MIN_CELLS ] || [ ! -z $GENES_MAX_CELLS ] ; then
    echo "-- filtering genes --"
    if [ -z $GENES_MIN_COUNTS ]; then GENES_MIN_COUNTS="0"; fi
    if [ -z $GENES_MAX_COUNTS ]; then GENES_MAX_COUNTS="0"; fi
    if [ -z $GENES_MIN_CELLS ]; then GENES_MIN_CELLS="0"; fi
    if [ -z $GENES_MAX_CELLS ]; then GENES_MAX_CELLS="0"; fi

    FULL_OUTPUT=$OUTPUT_BASENAME"_gene_filter.h5ad"
    python3 $SRC_PATH/filter_genes.py $DATA_FILE $FULL_OUTPUT $GENES_MIN_COUNTS $GENES_MAX_COUNTS $GENES_MIN_CELLS $GENES_MAX_CELLS
    exitOnError $? "Error during gene filtering."
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $NORMALIZE ] && [ "$NORMALIZE" -eq 1 ]; then
    echo "-- normalizing --"
    FULL_OUTPUT=$OUTPUT_BASENAME"_normalized.h5ad"
    python3 $SRC_PATH/generate_clusters_for_normalization.py $DATA_FILE
    exitOnError $? "Error generating clusters for normalization."
    Rscript $SRC_PATH/compute_size_factors.R
    exitOnError $? "Error computing size factors for normalization."
    python3 $SRC_PATH/normalize.py $DATA_FILE $FULL_OUTPUT
    exitOnError $? "Error normalizing the data."
    rm -f temp_clustered_for_scran.h5ad
    rm -f temp_size_factors.csv
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $BATCH_CORRECT ] && [ "$BATCH_CORRECT" -eq 1 ]; then
    echo "-- batch correcting --"
    FULL_OUTPUT=$OUTPUT_BASENAME"_batch_corrected.h5ad"
    echo "batch correction input: "$DATA_FILE
    echo "batch correction output: "$FULL_OUTPUT
    cp $DATA_FILE $FULL_OUTPUT
    Rscript $SRC_PATH/batch_correct.R $FULL_OUTPUT $BATCH_VAR
    exitOnError $? "Error during batch correction."
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $HIGH_VAR_GENES ]; then
    echo "-- selecting "$HIGH_VAR_GENES" high variance genes --"
    python3 $SRC_PATH/high_variance_genes.py $DATA_FILE $OUTPUT_BASENAME $HIGH_VAR_GENES
    exitOnError $? "Error filtering high variance genes."
    DATA_FILE=$OUTPUT_BASENAME"_high_variance_genes_subset.h5ad"
fi

if [ -z $COMPUTE_UMAP ]; then
    COMPUTE_UMAP=0
fi

if [ -z $COMPUTE_TSNE ]; then
    COMPUTE_TSNE=0
fi

if [[ "$COMPUTE_UMAP" -eq 1 || "$COMPUTE_TSNE" -eq 1 ]]; then
    FULL_OUTPUT=$OUTPUT_BASENAME"_dim_reduce.h5ad"
    python3 $SRC_PATH/dimension_reduction.py $DATA_FILE $FULL_OUTPUT $COMPUTE_UMAP $COMPUTE_TSNE
    exitOnError $? "Error creating umap or tsne plots."
    DATA_FILE=$FULL_OUTPUT
fi
