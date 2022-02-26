/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include"wallEvaCoupledMixed.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "rhoReactionThermo.H"

namespace Foam
{
namespace compressible
{
    // constructor from patch and internal field
    wallEvaCoupledMixedFvPatchScalarField::wallEvaCoupledMixedFvPatchScalarField
    (
        const fvPatch& p,
        const DimensionedField<scalar,volMesh>& iF
    ):
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    TnbrName_("undefined-Tnbr"),
    evaname_("undefined-evaname"),
    specieName_("undefined-specieName"),
    liquid_(nullptr),
    liquidDict_(nullptr),
    fluidside_(false),
    mass_(patch().size(),Zero),
    massOld_(patch().size(),Zero),
    mKdelta_(patch().size(),Zero),
    dmHfg_(patch().size(),Zero),
    Mcomp_(0.0),
    lastTimeStep_(0)
    {
        this->refValue() = 0.0;
        this->refGrad() = 0.0;
        this->valueFraction() = 1.0;
    }
    //constructor from patch ,internalfield and dictionary
    wallEvaCoupledMixedFvPatchScalarField::wallEvaCoupledMixedFvPatchScalarField
    (
        const fvPatch& p,
        const DimensionedField<scalar,volMesh>& iF,
        const dictionary& dict
    ):
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    Tname_(dict.lookupOrDefault<word>("Tnber","T")),
    evaname(dict.lookupOrDefault<word>("evaname","H2O")),
    specieName_(dict.lookupOrDefault<word>("specie","undefined-specieName")),
    fluidside_(false),
    liquid_(nullptr),
    liquidDict_(),
    mass_(patch().size(),Zero),
    massOld_(patch().size(),Zero),
    mKdelta_(patch().size(),Zero),
    dmHfg_(patch().size(),Zero),
    Mcomp_(dict.lookupOrDefault<scalar>("carrierMolWeight",0.0)),
    lastTimeStep_(0)
    {
        if (!isA<mappedPatchBase>(this->patch().patch()))
        {
            FatalErrorInFunction
                << "' not type '" << mappedPatchBase::typeName << "'"
                << "\n    for patch " << p.name()
                << " of field " << internalField().name()
                << " in file " << internalField().objectPath()
                << exit(FatalError);
        }

        if (dict.found("refValue"))
        {
            refValue() = scalarField("refValue", dict, p.size());
            refGrad() = scalarField("refGradient", dict, p.size());
            valueFraction() = scalarField("valueFraction", dict, p.size());
        }
        else
        {
            // Start from user entered data. Assume fixedValue.
            refValue() = *this;
            refGrad() = 0.0;
            valueFraction() = 1.0;
        }
        if (dict.found("specie"))
        {
            fluidside_ = true;
            Info << "   ============ found fluid side =========== "<<endl;
        }
        if (fluidside_)
        {
            liquidDict_ = dict.subDict("liquid");
            liquid_ = liquidProperties::New(liquidDict_.subDict(specieName_));
            lastTimeStep_ = patch().boundaryMesh().mesh().time().value();
        }
    }

    //map function autoMap
    void wallEvaCoupledMixedFvPatchScalarField::autoMap
    (
        const fvPatchFieldMapper& m
    )
    {
        mixedFvPatchField::autoMap(m)
        if (fluidside_)
        {
            mass_.autoMap(m);
            massOld_.autoMap(m);
            dmHfg_.autoMap(m);
        }        
    }

    //member function of update
    void wallEvaCoupledMixedFvPatchScalarField::updateCoeffs()
    {
        if (update())
        {
            return;
        }
        //tag
        int OldTag = UPstream::msgType();
        UPstream::msgType() = OldTag + 1;

        //get coupled information from the mappedPatchBase
        const mappedPatchBase & mpp = refCast<const mappedPatchBase>(patch().patch());
        //get the area of the patch face
        const scalarField & magSf = patch().magSf();
        //get the boundary mesh form patch
        const polyMesh & mesh = patch().boundaryMesh().mesh();
        //get the nbr mesh from mpp
        const polyMesh& nbrMesh = mpp.sampleMesh();
        //get the patchMesh index
        const label samplePatchi = mpp.sampleMesh().index();
        //get the nbrPatch from patch ID
        const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];
        //get the internal field using scalarField constructor
        scalarField TinternalField(patchInternalField());
        scalarField& Tpatch = *this;

        const volScalarField& TmeshField = static_cast<const volScalarField&>(internalField());
        scalarField TpatchOld(TmeshField.oldTime().boundaryField()[patch().index()]);

        //change typename
        typedef wallEvaCoupledMixedFvPatchScalarField thistype;
        //get nbr Tpatch from dict.TnbrName
        const fvPatchScalarField& nbrTp = nbrPatch.lookupPatchField<volScalarField,scalar>(TnberName_);
        //test check erro
        if (!isA<thisType>(nbrTp))
        {
            FatalErrorInFunction
                << "Patch field for " << internalField().name() << " on "
                << patch().name() << " is of type " << thisType::typeName
                << endl << "The neighbouring patch field " << TnbrName_ << " on "
                << nbrPatch.name() << " is required to be the same, but is "
                << "currently of type " << nbrTp.type() << exit(FatalError);
        }
        //change nbrTp into thistype
        const thistype& nbrField = refCast<const thistype&>(nbrTp);
        //get internal field from nbrpatch
        scalarField nbrInternalField(nbrField.patchInternalField());
        mpp.distribute(nbrInternalField);
        mKdelta_ = kappa(Tpatch)*patch().deltaCoeffs();
/*
        //does this is need?
        scalarField dm(patch().size(), Zero);
        scalarField dHspec(patch().size(), Zero);
        scalarField hPhaseChange(patch().size(), Zero);
        scalarField dmhPhaseChange(patch().size(), Zero);
        scalarField hRemovedMass(patch().size(), Zero);
        scalarField dmhRemoveMass(patch().size(), Zero);
*/      

        //get the mass changed 
        scalarField dm(patch().size(),Zero);
        //get the enthalpy in mass added
        scalarField hAddedMass(patch().size(),Zero);
        //get the unit mass changed enthalpy
        scalarField dmhAddedMass(patch().size(),Zero);

        if (fluidside_)
        {
            const scalar dt = mesh.time().deltaTValue();
            if (mesh.time().value() != lastTimeStep_)
            {
                //swap timestep's value from last time to next time
                lastTimeStep_ = mesh.time().value();
                massOld_ = mass_;
            }

            //get the thermo's porperties from the mesh
            const rhoReactionThermo& thermo = mesh.lookupObject<rhoReactionThermo>(basicThermo::dictName);
            //get the composition of the thermo model
            const basicSpecieMixture& composition = thermo.composition();
            //get the species index of the composition
            label specieIndex = composition.species()[specieName_];

            //set the mass change temp
            volScalarField& MassSource = const_cast<volScalarField&>(mesh.lookupObject<volScalarField>(specieName_));
            //set the energy change temp
            volScalarField& EnergySource = const_cast<volScalarField&>(mesh.lookupObject<volScalarField>(specieName_));

            //get delta=1/distance cell center to patch
            const scalarField myDelta_(patch().deltaCoeffs());
            //get the cp 
            scalarField cp(patch().size(),Zero);
            //get the rho
            scalarField liquidRho(patch().size(),Zero);
            //get the mass species of the patch
            const fvPatchScalarField& Ypatch = patch().lookupPatchField<volScalarField,scalar>(specieName_);
            //get the pressure of the patch
            const fvPatchScalarField& pPatch = patch().lookupPatchField<volScalarField,scalar>("p");
            //get the denisty of the patch
            const fvPatchScalarField& rhoPatch = patch().lookupPatchField<volScalarField,scalar>("rho");
            //get the turbulent viscosity of the patch 
            const fvPatchScalarField& nutPatch = patch().lookupPatchField<volScalarField,scalar>("nut");
            //get the liquid viscosity of the patch (maybe not)
            //const fvPatchScalarField& muPatch = patch().lookupPatchField<volScalarField,scalar>("thermo:mu");
            //internal field get Y p rho in cell center *** not the patch
            const scalarField Yinternal(Ypatch.patchInternalField());
            const scalarField pInternal(pPatch.patchInternalField());
            const scalarField rhoInternal(rhoPatch.patchInternalField());
            //get the face cell's index list
            const labelList& faceCells = patch().faceCells();
            //because this env is standred gas-solid ,so the Sc = 0.7
            const scalar Sct = 0.7;
            const scalar Sc = 0.7;

            //caculate the boundary cell value
            forAll(Tpatch,faceI)
            {
                const scalar Tface = Tpatch[faceI];
                const scalar Tcell = Tinternal[faceI];
                const scalar pFace = pPatch[faceI];
                const scalar pcell = pInternal[faceI];

                //const scalar muFace = muPatch[faceI];
                //set the liquid viscosity in nuCell is 2.95e-5 (kinematic viscosity)?
                const scalar nuCell = 2.95e-5;
                //set the rho of face
                const scalar rhoFace = rhoPatch[faceI];
                //set the rho of cell
                const scalar rhoCell = rhoInternal[faceI];
                //set the liquid (Dynamic viscosity)?
                const scalar muCell = nuCell * rhoCell;
                //set the turbulent liquid (Dynamic viscosity)?
                const scalar mutFace = nutPatch[faceI]*rhoFace;
                //get the saturation of the cell
                const scalar pSatCell = liquid_->pv(pFace, Tcell);
                //get the saturation of the face
                const scalar pSatFace = liquid_->pv(pFace, Tface);
                //get the molecular | Molar weight [Kg/Kmol]
                const scalar Mv = liquid_->W();
                //get the species of the cell
                const scalar Ycell = Yinternal[faceI];
                //get the 1/distance of the face
                const scalar deltaFace = myDelta[faceI];

                //get the cpacity of the liquid 
                cp[faceI] = liquid_->Cp(pFace,Tface);
                //get the phase change energy
                hPhaseChange[faceI] = liquid_->hl(pFace,Tface);
                //get the mass energy from phase
                hRemovedMass[faceI] = composition.Hs(specieIndex,pcell,Tcell);


                const scalar YsatFace = pSatFace/pFace*Mv/Mcomp_;
                const scalar gamma = muCell + mutFace / Sct;
                //calculate mass add
                dm[faceI] = gamma * deltaFace * ( Ycell - YsatFace )/(1 - YsatFace);
                //calculate the liquid rho
                liquidRho[faceI] = liquid_->rho(pFace,Tface);                
            }

            scalarField &massFluxOut = outputScalarField(specieName_ + "MassFlux",dimMass/dimArea/dimTime,refCast<const fvMesh>(mesh)).boundaryFieldRef()[patch().index()];
            //define a new avarible to store massflux
            massFluxOut = dm;

            //avariable for liquidRho
            scalarField &rhoLiquid = outputScalarField("densityLiquid",dimDensity,refCast<const fvMesh>(mesh));
            rhoLiquid = rhoPatch.internalField()*Ypatch.internalField();

            //avariable for gasRho
            scalarField &rhoGas = outputScalarField("densityGas",dimDensity,refCast<const fvMesh>(mesh));

            //get phase of air from dictfile
            const volScalarField& Yair = mesh.lookupObject<volScalarField>("AIR");

            rhoGas = rhoPatch.internalField()*Yair;
            //calculate mass change of phase; and compare with 0,makesure +
            mass_ = massOld_ * dm * dt * magSf;
            mass_ = max(mass_,scalar(0));
            //calculate heat flux change of phase
            dmHfg_ = dm * hPhaseChange;
            dHspec = dm * hRemovedMass;

            forAll(faceCells,faceI)
            {
                const label cellI = faceCells[faceI];
                MassSource[cellI] = -dm[faceI] * magSf[faceI] / mesh.cellVolumes()[faceI];
                EnergySource[cellI] = -dm[faceI] * magSf[faceI] / mesh.cellVolumes()[faceI]*hRemovedMass[faceI];
            }
        }//end liquidside calculate

        //calclulate k*delta by resistence
        scalarField KDeltaNbr;





        mixedFvPatchScalarField::updateCoeffs();

        if (fluidside_)
        {
            scalar Qdm = gSum(dm*magSf)
        }
        
        
    }


}
}






