// Initialise solid field pointer lists
PtrList<coordinateSystem> coordinates(solidRegions.size());
PtrList<solidThermo> thermos(solidRegions.size());
PtrList<radiation::radiationModel> radiations(solidRegions.size());
PtrList<fv::options> solidHeatSources(solidRegions.size());
PtrList<volScalarField> betavSolid(solidRegions.size());
PtrList<volSymmTensorField> aniAlphas(solidRegions.size());

//add to solid section
PtrList<volScalarField> interSolid(solidRegions.size());
PtrList<volScalarField> cellVolume(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info << "   Adding to thermoSolid\n"<<endl;
    thermoSolid.set(i,heThermo(fluidRegions[i],world::null));

    Info<< "    Adding to thermos\n" << endl;
    thermos.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to radiations\n" << endl;
    radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

    //add the internal 
    Info<< "    Addomg the interface \n"<<endl;
    interSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "H2O",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("H2O",dimless,0)
        )
    )
    //add the cell volume to the file
    Info<<"     Adding the cell volume \n"<<endl;
    cellVolume.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "cellVolume",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
            solidRegions[i],
            dimensionedScalar("temp",dimless,0)
        )
    )
    forAll(cellVolume[i],cellI)
    {
        cellVolume[i][cellI] = solidRegions[i].V()[cellI];//add the cell volume to the cellvolume file 
    }


    if (!thermos[i].isotropic())
    {
        Info<< "    Adding coordinateSystems\n" << endl;
        coordinates.set
        (
            i,
            coordinateSystem::New
            (
                solidRegions[i],
                thermos[i],
                coordinateSystem::typeName_()
            )
        );

        tmp<volVectorField> tkappaByCp =
            thermos[i].Kappa()/thermos[i].Cp();

        aniAlphas.set
        (
            i,
            new volSymmTensorField
            (
                IOobject
                (
                    "Anialpha",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedSymmTensor(tkappaByCp().dimensions(), Zero),
                zeroGradientFvPatchSymmTensorField::typeName
            )
        );

        aniAlphas[i].primitiveFieldRef() =
            coordinates[i].transformPrincipal
            (
                solidRegions[i].cellCentres(),
                tkappaByCp()
            );
        aniAlphas[i].correctBoundaryConditions();

    }

    IOobject betavSolidIO
    (
        "betavSolid",
        runTime.timeName(),
        solidRegions[i],
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (betavSolidIO.typeHeaderOk<volScalarField>(true))
    {
        betavSolid.set
        (
            i,
            new volScalarField(betavSolidIO, solidRegions[i])
        );
    }
    else
    {
        betavSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "betavSolid",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar("1", dimless, scalar(1))
            )
        );
    }
}