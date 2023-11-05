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

auto hullA = hull_a[i];
auto hullB = hull_b[i];
auto branA = bran_a[i];
auto branB = bran_b[i];
auto endoA = endo_a[i];
auto endoB = endo_b[i];

forAll(mesh.cellZones(),cellZone)
{
    const labelList& selectedCells(mesh.cellZones()[cellZone]);
    const word ZoneName(mesh.cellZones()[cellZone].name());
    if (ZoneName == "hull")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = hullA*Foam::exp(-hullB/solidTmp[cellIndex])/3600.0;
        }
    }else if (ZoneName == "bran")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = branA*Foam::exp(-branB/solidTmp[cellIndex])/3600.0;
        }
    }else if (ZoneName == "endosperm")
    {
        forAll(selectedCells, loopIndex)
        {
            const label& cellIndex = selectedCells[loopIndex];
            diffusionSolid[i][cellIndex] = endoA*Foam::exp(-endoB/solidTmp[cellIndex])/3600.0;
        }
    }
}

diffusionSolid[i].boundaryFieldRef() = 484.0/3600.0*Foam::exp(-7380.8/thermos[i].T().boundaryFieldRef());
diffusionSolid[i].correctBoundaryConditions();

//add
volScalarField& interF = interSolid[i];
volScalarField& Dsolid = diffusionSolid[i];
volScalarField& kappa = solidKappa[i];
volScalarField& cp_g = Cp_g[i];
fv::options& fvOptions = solidHeatSources[i];
//calculated the gradient of moisture and temperature
YiGradient[i] = fvc::grad(interF);
TemGradient[i] = fvc::grad(solidTmp);
