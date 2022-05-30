fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];

volScalarField& Tsolidtmp = thermo.T();

volScalarField& interEFF_coe = interEffSolid[i];

{
    //for loop internal cell value
    forAll(mesh.cellZones(),cellZone)
    {
        const labelList & selectedCells(mesh.cellZones()[cellZone]);
        forAll(selectedCells,loopindex)
        {
            const label& cellIndex = selectedCells[loopindex];
            diffusionSolid[i][cellIndex] = interEFF_coe[cellIndex] * Foam::exp(-27.36 / (0.008314*Tsolidtmp[cellIndex]));
        }
    }
    
    //for loop boundary cell value
    volScalarField::Boundary& Tbf = Tsolidtmp.boundaryFieldRef();
    volScalarField::Boundary& dbf = diffusionSolid[i].boundaryFieldRef();
    volScalarField::Boundary& ibf = interEFF_coe.boundaryFieldRef();

    forAll(mesh.boundary(),patch)
    {
        const word& patchName = mesh.boundary()[patch].name();
        Info << " patch's name:"<<patchName << endl;
        forAll(mesh.boundary()[patch],faceI)
        {
            dbf[patch][faceI] = ibf[patch][faceI] * Foam::exp(-27.36 / (0.008314*Tbf[patch][faceI]));
        }
    }
}
