#!/bin/sh

cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

mpirun -np 16 kernelPhaseTrasition4FOAM -parallel &> run.log &
