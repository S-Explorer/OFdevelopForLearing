#!/bin/sh

cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

rm -rf log.*
rm -rf constant/*/polyMesh

foamCleanPolyMesh
echo "clean mesh done !"

blockMesh > log.blockMesh
echo "blockMesh done !"

snappyHexMesh -overwrite > log.snappy
echo "snappyHexMesh done !"

splitMeshRegions -cellZones -overwrite > log.splitMeshRegions
echo "splitMeshRegions done !"

topoSet -region hull > log.topoSet
topoSet -region fluid >> log.ropoSet
echo "topoSet cell Zones done !"

checkMesh -allRegions > log.checkMesh
echo "checkMesh"

decomposePar -allRegions > log.decompose
echo "decomposePar done !"