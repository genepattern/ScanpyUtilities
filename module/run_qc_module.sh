#!/bin/bash

# version information
echo "Version v0.1"

# default values for parameters
PY_EXEC="python"
SRC_PATH=.

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

    esac
done

if [ ! -z ${GENOME} ] && [ ${DATA_FILE: -3} == ".h5" ]; then
    echo "-- converting from H5 to H5AD format --"
    FULL_OUTPUT="${OUTPUT_BASENAME}.h5ad"
    eval "${PY_EXEC} ${SRC_PATH}/convert_to_h5ad.py ${DATA_FILE} ${GENOME} ${FULL_OUTPUT}"
    exitOnError $? "Error during conversion from h5 to h5ad format."
    DATA_FILE=${FULL_OUTPUT}
fi

if [ -z ${MITO_FILE} ]; then MITO_FILE="SKIP"; fi

echo "-- adding default annotations --"
eval "${PY_EXEC} ${SRC_PATH}/qc_experiment.py ${DATA_FILE} ${MITO_FILE} ${OUTPUT_BASENAME}"
exitOnError $? "Error during quality control analysis."
