fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];
Info<<" \t *** set SolidRegions rho filed !" << endl;
tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

Info<<" \t *** set SolidRegions cp filed !" << endl;
tmp<volScalarField> tcp = thermo.Cp();
const volScalarField& cp = tcp();

tmp<volSymmTensorField> taniAlpha;
if (!thermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCp = thermo.Kappa()/cp;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.transformPrincipal
        (
            mesh.cellCentres(),
            tkappaByCp()
        );

    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}

Info<<" \t *** set SolidRegions h filed !" << endl;
volScalarField& h = thermo.he();
//add
Info<<" \t *** set SolidRegions specis filed !" << endl;
volScalarField& interF = interSolid[i];

Info<<" \t *** set SolidRegions source filed !" << endl;
volScalarField& HSource = sourceHeat[i];

Info<<" \t *** set SolidRegions diffusion filed !" << endl;
volScalarField& Dsolid = diffusionSolid[i];

const volScalarField& betav = betavSolid[i];

Info<<" \t *** set SolidRegions fvOptions !" << endl;
fv::options& fvOptions = solidHeatSources[i];

Info<<"-------------end set SolidRegions !" << endl;
