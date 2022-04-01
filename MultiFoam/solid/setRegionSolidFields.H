fvMesh& mesh = solidRegions[i];
solidThermo& thermo = thermos[i];

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

volScalarField HSource
(
		IOobject
	(
			"HSource",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
	),
	mesh,
	dimensionedScalar("HSource",dimPower/dimVolume,Zero)
);

volScalarField Dsolid
(
    IOobject
    (
            "Dsolid",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimViscosity,Zero)
);

if (!Dsolid.headerOk())
{    
    IOdictionary solidtransProperties
    (
        IOobject
        (
            "SolidProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    Dsolid = dimensionedScalar("Dsolid",dimViscosity,solidtransProperties);
}

volScalarField& h = thermo.he();
//add
volScalarField& interF = interSolid[i];

const volScalarField& betav = betavSolid[i];

fv::options& fvOptions = solidHeatSources[i];
