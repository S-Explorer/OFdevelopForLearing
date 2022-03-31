/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenCFD Ltd.
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

#ifndef massTransferCoupledFvPatchScalarField_H
#define massTransferCoupledFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "liquidProperties.H"
#include "autoPtr.H"

namespace Foam
{
namespace compressible{
/*---------------------------------------------------------------------------*\
           Class massTransferCoupledFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class massTransferCoupledFvPatchScalarField:public mixedFvPatchScalarField,public temperatureCoupledBase
{
private:
    //蒸发速率计算
    scalar EvaporationRate(scalar) const;
    //导出物理场
    volScalarField& outputScalarField
    (
        const word& fieldName,
        const dimensionSet& dimSet,
        const fvMesh& mesh
    );
    //压力场
    const word pName_;

    //速度场
    const word UName_;

    //密度场
    const word rhoName_;

    //动力粘度
    const word muName_;

    //邻区的温度场
    const word TnbrName_;

    //传质的物质（default ： h2o）
    const word specieName_;

    //能量方程中的源项
    const word sourceTerm_;

    //特征长度
    scalar L_;

    //舍伍德数
    scalar Sh_;

    //获取固体侧的扩散系数
    const word SolidDivCoe_;

    //流体属性
    autoPtr<liquidProperties> liquid_;
    dictionary liquidDict_;

    //面的质量积累
    scalarField mass_;
    scalarField massOld_;

    //汽化温度
    scalar Tvap_;

    //kappa * deltacoeffs
    scalarField myKDelta_;

    //相变能
    scalarField dmHfg_;

    //热惯性
    scalarField mpCpTp_;

    //平均分子质量-载体空气
    scalar Mcomp_;

    //计算侧（流体 or 固体）
    bool fluid_;

    //比热容
    scalarField cp_;

    //密度
    scalarField rho_;

    //- Last time step for which mass of fluid film was calculated
    scalar lastTimeStep_;

public:
    //public functions
    TypeName("massTransferCoupled");

    //构造函数
    massTransferCoupledFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>& 
    );

    massTransferCoupledFvPatchScalarField
    (
        const fvPatch& ,
        const DimensionedField<scalar, volMesh>& ,
        const dictionary&
    );
    
    massTransferCoupledFvPatchScalarField
    (
        const massTransferCoupledFvPatchScalarField&,
        const fvPatch& ,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    massTransferCoupledFvPatchScalarField
    (
        const massTransferCoupledFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    virtual void updateCoeffs();

    virtual void write(Ostream&) const;

};



}
}
#endif
