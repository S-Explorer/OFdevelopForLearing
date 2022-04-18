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
PtrList<volScalarField> sourceHeat(solidRegions.size());
PtrList<volScalarField> diffusionSolid(solidRegions.size());
PtrList<IOdictionary> dictionSolid(solidRegions.size());

// Populate solid field pointer list
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermos\n" << endl;
    thermos.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to radiations\n" << endl;
    radiations.set(i, radiation::radiationModel::New(thermos[i].T()));

    Info<< "    Adding to fvOptions\n" <<endl;
    solidHeatSources.set
    (
        i,
        new fv::options(solidRegions[i])
    );

    //add the internal 
    Info<< "    Adding the interface \n"<<endl;
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
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("H2O",dimless,0)
        )
    );

    //setting for heat source 
    Info << "   Adding the heat source change tmp\n" << endl;
    sourceHeat.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "HSource",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],  
            dimensionedScalar("HSource",dimPower/dimVolume,Zero)          
        )
    );

    //read the specis diffusion number
    Info << "   Adding the diffusionSolid number \n" << endl;
    diffusionSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Dsolid",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            solidRegions[i],            
            dimensionedScalar("Dsolid",dimViscosity,Zero)
        )
    );
    //add the properties file in solid region constant folder!
    dictionSolid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "SolidProperties",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    );

    diffusionSolid[i] = dimensionedScalar("Dsolid",dimViscosity,dictionSolid[i]);

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
            ),
            solidRegions[i],
            dimensionedScalar("cellVolume",dimVolume,0)
        )
    );

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
