// Initialise solid field pointer lists
PtrList<solidThermo> thermos(solidRegions.size());
PtrList<radiation::radiationModel> radiations(solidRegions.size());
PtrList<fv::options> solidHeatSources(solidRegions.size());
//add to solid section
PtrList<volScalarField> interSolid(solidRegions.size());
PtrList<volScalarField> cellVolume(solidRegions.size());
PtrList<volScalarField> sourceHeat(solidRegions.size());
PtrList<volScalarField> diffusionSolid(solidRegions.size());
PtrList<IOdictionary> dictionSolid(solidRegions.size());
// solid thermo physical properties
PtrList<volScalarField> solidCp(solidRegions.size());
PtrList<volScalarField> solidKappa(solidRegions.size());
PtrList<volVectorField> YiGradient(solidRegions.size());
PtrList<volVectorField> TemGradient(solidRegions.size());

// Populate solid field pointer lists
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
                IOobject::AUTO_WRITE
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

    dimensionedScalar TemOnes("ones",thermos[i].T().dimensions(),1.0);
    tmp<volScalarField> tmpdiffusion = 3.01e7*Foam::exp(-3510*TemOnes/thermos[i].T());

    diffusionSolid[i].primitiveFieldRef() = tmpdiffusion();
    diffusionSolid[i].boundaryFieldRef() = 3.01e7*Foam::exp(-3510/thermos[i].T().boundaryFieldRef());
    diffusionSolid[i].correctBoundaryConditions();
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

    YiGradient.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "H2OGradient",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedVector("zero",dimless/dimLength,vector::zero)
        )
    );

    YiGradient[i] = fvc::grad(interSolid[i]);

    TemGradient.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "TemperatureGradient",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedVector("zero",dimTemperature/dimLength,vector::zero)
        )
    );

    TemGradient[i] = fvc::grad(thermos[i].T());

    Info << "   Adding the solid heat capacity \n" << endl;
    solidCp.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "solidCp",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("zero",dimEnergy/dimMass/dimTemperature,Zero)
        )
    );

    Info << "   Adding the solid kappa \n" << endl;
    solidKappa.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "solidKappa",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("zero",dimEnergy/dimLength/dimTime/dimTemperature,Zero)
        )
    );

    Info << "update the solid kappa ! \n" << endl;

    tmp<volScalarField> tmpKappa = (0.0637+0.0958*interSolid[i])/(0.656-0.475*interSolid[i]);

    solidKappa[i].primitiveFieldRef() = tmpKappa();
    solidKappa[i].boundaryFieldRef() = (0.0637+0.0958*interSolid[i].boundaryFieldRef())/(0.656-0.475*interSolid[i].boundaryFieldRef());
    solidKappa[i].correctBoundaryConditions();

    tmp<volScalarField> tmpCp = thermos[i].Cp();
    solidCp[i] = tmpCp();
}
