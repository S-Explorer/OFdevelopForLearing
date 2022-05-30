fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();
tmp<volScalarField> tcp = thermo.Cp();
const volScalarField& cp = tcp();

// calculate the solid kappa value 
tmp<volScalarField> tmpKappa = (0.0637+0.0958*interSolid[i])/(0.656-0.475*interSolid[i]);
solidKappa[i].primitiveFieldRef() = tmpKappa();
solidKappa[i].boundaryFieldRef() = (0.0637+0.0958*interSolid[i].boundaryFieldRef())/(0.656-0.475*interSolid[i].boundaryFieldRef());
solidKappa[i].correctBoundaryConditions();
//get ther solidthermo 's enathlty and temperature
volScalarField& h = thermo.he();
volScalarField& solidTmp = thermo.T();
//calculate ther solid diffusion value 
dimensionedScalar TemOnes("ones",thermos[i].T().dimensions(),1.0);

forAll(mesh.cellZones(),cellZone)
{
    const labelList& selectedCells(mesh.cellZones()[cellZone]);
    const word ZoneName(mesh.cellZones()[cellZone].name());
    if (ZoneName == "hull")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = max(3.07e-7*Foam::exp(-3510/solidTmp[cellIndex]),VSMALL);
        }
    }else if (ZoneName == "bran")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = max(1.77e-4*Foam::exp(-5390/solidTmp[cellIndex]),VSMALL);            
        }
    }else if (ZoneName == "endosperm")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = max(2.58e-3*Foam::exp(-5680/solidTmp[cellIndex]),VSMALL);            
        }
    }
}

diffusionSolid[i].boundaryFieldRef() = 3.07e-7*Foam::exp(-3510/thermos[i].T().boundaryFieldRef());
diffusionSolid[i].correctBoundaryConditions();

//add
volScalarField& interF = interSolid[i];
volScalarField& HSource = sourceHeat[i];
volScalarField& Dsolid = diffusionSolid[i];
volScalarField& kappa = solidKappa[i];
fv::options& fvOptions = solidHeatSources[i];
//calculated the gradient of moisture and temperature
YiGradient[i] = fvc::grad(interF);
TemGradient[i] = fvc::grad(solidTmp);
