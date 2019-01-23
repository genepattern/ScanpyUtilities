#!/bin/bash

for i in "$@"
do
    case $i in
        --data.file=*)
        DATA_FILE="${i#*=}"
        ;;

        --output.basename=*)
        OUTPUT_BASENAME="${i#*=}"
        ;;

        --annotate=*)
        ANNOTATE="${i#*=}"
        ;;

        --cells.min.counts=*)
        CELLS_MIN_COUNTS="${i#*=}"
        ;;

        --cells.max.counts=*)
        CELLS_MAX_COUNTS="${i#*=}"
        ;;

        --cells.min.genes=*)
        CELLS_MIN_GENES="${i#*=}"
        ;;

        --cells.max.genes=*)
        CELLS_MAX_GENES="${i#*=}"
        ;;

        --genes.min.counts=*)
        GENES_MIN_COUNTS="${i#*=}"
        ;;

        --genes.max.counts=*)
        GENES_MAX_COUNTS="${i#*=}"
        ;;

        --genes.min.cells=*)
        GENES_MIN_CELLS="${i#*=}"
        ;;

        --genes.max.cells=*)
        GENES_MAX_CELLS="${i#*=}"
        ;;

        --normalize=*)
        NORMALIZE="${i#*=}"
        ;;

        --batch.correct=*)
        BATCH_CORRECT="${i#*=}"
        ;;

        --batch.variable=*)
        BATCH_VAR="${i#*=}"
        ;;

        --n.high.variance.genes=*)
        HIGH_VAR_GENES="${i#*=}"
        ;;

        *)
            echo "ERROR: unrecognized option"
            echo "${i%=*}"
            exit 1
        ;;
    esac
done

if [ ! -z $ANNOTATE ] && [ "$ANNOTATE" -eq 1 ]; then
    echo "-- adding default annotations --"
    FULL_OUTPUT=$OUTPUT_BASENAME"_annotated.h5ad"
    python3 add_default_annotations.py $DATA_FILE $FULL_OUTPUT
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $CELLS_MIN_COUNTS ] || [ ! -z $CELLS_MAX_COUNTS ] || [ ! -z $CELLS_MIN_GENES ] || [ ! -z $CELLS_MAX_GENES ] ; then
    echo "-- filtering cells --"
    if [ -z $CELLS_MIN_COUNTS ]; then CELLS_MIN_COUNTS="0"; fi
    if [ -z $CELLS_MAX_COUNTS ]; then CELLS_MAX_COUNTS="0"; fi
    if [ -z $CELLS_MIN_GENES ]; then CELLS_MIN_GENES="0"; fi
    if [ -z $CELLS_MAX_GENES ]; then CELLS_MAX_GENES="0"; fi

    FULL_OUTPUT=$OUTPUT_BASENAME"_cell_filter.h5ad"
    python3 filter_cells.py $DATA_FILE $FULL_OUTPUT $CELLS_MIN_COUNTS $CELLS_MAX_COUNTS $CELLS_MIN_GENES $CELLS_MAX_GENES
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $GENES_MIN_COUNTS ] || [ ! -z $GENES_MAX_COUNTS ] || [ ! -z $GENES_MIN_CELLS ] || [ ! -z $GENES_MAX_CELLS ] ; then
    echo "-- filtering genes --"
    if [ -z $GENES_MIN_COUNTS ]; then GENES_MIN_COUNTS="0"; fi
    if [ -z $GENES_MAX_COUNTS ]; then GENES_MAX_COUNTS="0"; fi
    if [ -z $GENES_MIN_CELLS ]; then GENES_MIN_CELLS="0"; fi
    if [ -z $GENES_MAX_CELLS ]; then GENES_MAX_CELLS="0"; fi

    FULL_OUTPUT=$OUTPUT_BASENAME"_gene_filter.h5ad"
    python3 filter_genes.py $DATA_FILE $FULL_OUTPUT $GENES_MIN_COUNTS $GENES_MAX_COUNTS $GENES_MIN_CELLS $GENES_MAX_CELLS
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $NORMALIZE ] && [ "$NORMALIZE" -eq 1 ]; then
    echo "-- normalizing --"
    FULL_OUTPUT=$OUTPUT_BASENAME"_normalized.h5ad"
    python3 generate_clusters_for_normalization.py $DATA_FILE
    Rscript compute_size_factors.R
    python3 normalize.py $DATA_FILE $FULL_OUTPUT
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $BATCH_CORRECT ] && [ "$BATCH_CORRECT" -eq 1 ]; then
    echo "-- batch correcting --"
    FULL_OUTPUT=$OUTPUT_BASENAME"_batch_corrected.h5ad"
    echo "batch correction input: "$DATA_FILE
    echo "batch correction output: "$FULL_OUTPUT
    cp $DATA_FILE $FULL_OUTPUT
    Rscript batch_correct.R $FULL_OUTPUT $BATCH_VAR
    DATA_FILE=$FULL_OUTPUT
fi

if [ ! -z $HIGH_VAR_GENES ]; then
    echo "-- selecting "$HIGH_VAR_GENES" high variance genes --"
    python3 high_variance_genes.py $DATA_FILE $OUTPUT_BASENAME $HIGH_VAR_GENES
fi
