echo "clean Old Mesh !"
rm -rf log.*
rm -rf constant/fluid/polyMesh constant/kernel/polyMesh

foamCleanPolyMesh
echo "clean mesh done !"

blockMesh |tee log.blockMesh
echo "blockMesh done ! > log.blockMesh"

snappyHexMesh -overwrite | tee log.snappy
echo "snappyHexMesh done ! > log.snappy"

rm -rf 0/cell* 0/thick* 0/point* 0/nSurfaceLayers

echo "generate 2D mesh in this region !"
extrudeMesh -overwrite |tee log.extrudeMesh
echo "extrudeMesh done ! > log.extrude"

splitMeshRegions -cellZones -overwrite | tee log.splitMesh
echo "splitMesh Done ! > log.splitMesh" 

checkMesh -allRegions | tee log.checkMesh

for i in fluid kernel
do
	cp constant/${i}/boundary constant/${i}/polyMesh
done

rm -rf 1***
echo "================== boundary define done !==================="