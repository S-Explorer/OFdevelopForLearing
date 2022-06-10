tmpKappa = (0.0637+0.0958*interSolid[i])/(0.656-0.475*interSolid[i]);
solidKappa[i].primitiveFieldRef() = tmpKappa();
solidKappa[i].boundaryFieldRef() = (0.0637+0.0958*interSolid[i].boundaryFieldRef())/(0.656-0.475*interSolid[i].boundaryFieldRef());
solidKappa[i].correctBoundaryConditions();

forAll(mesh.cellZones(),cellZone)
{
    const labelList& selectedCells(mesh.cellZones()[cellZone]);
    const word ZoneName(mesh.cellZones()[cellZone].name());
    if (ZoneName == "hull")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = 3.07e-7*Foam::exp(-3510/solidTmp[cellIndex]);
        }
    }else if (ZoneName == "bran")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = 1.77e-4*Foam::exp(-5390/solidTmp[cellIndex]);            
        }
    }else if (ZoneName == "endosperm")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = 2.58e-3*Foam::exp(-5680/solidTmp[cellIndex]);            
        }
    }
}

diffusionSolid[i].boundaryFieldRef() = 3.07e-7*Foam::exp(-3510/thermos[i].T().boundaryFieldRef());
diffusionSolid[i].correctBoundaryConditions();