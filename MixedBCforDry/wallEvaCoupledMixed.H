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

useage:
    this boundary condition is heatAndMass transfer in solid-gas system,
    which deal a simple trouble about vapor acrose the interface between solid and gas;

    patchname
    {
        type    compressible::wallEvaCoupledMixed;
        Tnbr    T;
        evaname H2O;
        specie  H2O;

    }

\*---------------------------------------------------------------------------*/


#ifndef _WALLEVACOUPLEDMIXED_
#define _WALLEVACOUPLEDMIXED_

#include "mixedFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "scalarList.H"
#include "volFields.H"
#include "liquidProperties.H"

namespace Foam
{
namespace compressible
{
/*---------------------------------------------------------------------------*\
                Class wallEvaCoupledMixed Declaration
\*---------------------------------------------------------------------------*/
class wallEvaCoupledMixedFvPatchScalarField:
public mixedFvPatchScalarField,
public temperatureCoupledBase
{
private:
    //name of field on the neighbour region
    const word Tnbrname_;
    //name of the species on which the mass transferred (default H2O)
    const word specieName_;
    //name of the H2O evap filed's name 
    const word evaname_;
    //liquid properties
    autoPtr<liquidProperties> liquid_;
    //liquid dictionary
    dictionary liquidDict_;
    //mass accumulated on face 
    scalarField mass_;
    //mass accumulated on face in last time step
    scalarField massOld_;

    //fluid side ?
    bool fluidside_;
    //K
    scalarField mKdelta_;
    //phase change energy
    scalarField dmHfg_;
    //last time step for which mass of fluid film
    scalar lastTimeStep_;

public:
    //type infomtion
    TypeName("compressible::wallEvaCoupledMixed");
    //constructors
    // 1 constructors from patch and internal field
    wallEvaCoupledMixedFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar,volMesh>&
    );
    // 2 constructors from patch ,internal field and dictionary
    wallEvaCoupledMixedFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar,volMesh>&,
        const dictionary&
    );
    // 3 constructors from patch //delete


    //member function
    //map from self given a mapping object
    virtual void autoMap(const fvPatchFieldMapper&);
    //reverse

    //return selfK

    //return dmHfg
    const scalarField dmHfg() const
    {
        return dmHfg_;
    }

    //update the coefficient
    virtual void updateCoeffs();

    //write the reault file 
    virtual void write(Ostream&) const;

}
}
}
#endif