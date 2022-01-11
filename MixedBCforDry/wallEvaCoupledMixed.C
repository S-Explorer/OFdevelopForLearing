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
            const scalar Sc = 0.7;

            //caculate the boundary cell value
            forAll(Tpatch,faceI)
            {
                const scalar Tface = Tpatch[faceI];
                const scalar Tcell = Tinternal[faceI];
                const scalar pFace = pPatch[faceI];
                const scalar pcell = pInternal[faceI];
                
            }


            
        }




        mixedFvPatchScalarField::updateCoeffs();
        
    }


}
}






