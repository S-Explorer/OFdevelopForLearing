fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];

CompositionModel<heThermo>& mixture = mixSolid[i];
heThermo& thermo = mixture.thermo();
basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y(); 

tmp<volScalarField> trho = thermo.rho();
const volScalarField& rho = trho();

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


volScalarField& h = thermo.he();
//add
volScalarField& interF = interSolid[i];

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];
