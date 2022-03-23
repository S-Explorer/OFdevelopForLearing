/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "massTransferCoupledFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mappedPatchBase.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam{
namespace compressible{
//计算蒸发率
scalar massTransferCoupledFvPatchScalarField::EvaporationRate(scalar T) const
{
    if(T < scalar(430))
    {
        scalar A = 5.13e10;//Amoi   /s
        scalar E = 88000;//Emoi     J/mol
        scalar R = 8.31446261815;//R     J/kmol
        scalar k = A * exp(-1 * E / R / T);
        return k; 
    }else{
        return scalar(0.9);
    }
    return scalar(0);    
}
//导出物理场
volScalarField& massTransferCoupledFvPatchScalarField::outputScalarField
(
    const word& fieldName,
    const dimensionSet& dimSet,
    const fvMesh& mesh
)
{
    if (!mesh.foundObject<volScalarField>(fieldName))
    {
        tmp<volScalarField> tField
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimSet, Zero)
            )
        );

        tField.ptr()->store();
    }

    return const_cast<volScalarField &>(mesh.lookupObject<volScalarField>(fieldName));
}

massTransferCoupledFvPatchScalarField::massTransferCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
):
mixedFvPatchScalarField(p,iF),
temperatureCoupledBase(patch(),"undefined","undefined","undefined-K","undefined-alpha"),
pName_("p"),
UName_("U"),
rhoName_("rho"),
muName_("thermo:mu"),
TnbrName_("T"),
specieName_("none"),
sourceTerm_("eva"),
SolidDivCoe_("Ds"),
liquid_(nullptr),
liquidDict_(nullptr),
mass_(patch().size(), Zero),
massOld_(patch().size(), Zero),
Tvap_(0.0),
myKDelta_(patch().size(), Zero),
dmHfg_(patch().size(), Zero),
mpCpTp_(patch().size(), Zero),
Mcomp_(0.0),
fluid_(false),
cp_(patch().size(), Zero),
rho_(patch().size(), Zero),
lastTimeStep_(0.0)
{
    this->refValue() = 0.0 ;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}

massTransferCoupledFvPatchScalarField::massTransferCoupledFvPatchScalarField
(
    const massTransferCoupledFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
):
mixedFvPatchScalarField(psf,p,iF,mapper),
temperatureCoupledBase(patch(),psf),
pName_(psf.pName_),
UName_(psf.UName_),
rhoName_(psf.rhoName_),
muName_(psf.muName_),
TnbrName_(psf.TnbrName_),
specieName_(psf.specieName_),
sourceTerm_(psf.sourceTerm_),
SolidDivCoe_(psf.SolidDivCoe_),
liquid_(psf.liquid_.clone()),
liquidDict_(psf.liquidDict_),
mass_(psf.mass_,mapper),
massOld_(psf.massOld_,mapper),
Tvap_(psf.Tvap_),
myKDelta_(psf.myKDelta_,mapper),
dmHfg_(psf.dmHfg_,mapper),
mpCpTp_(psf.mpCpTp_,mapper),
Mcomp_(psf.Mcomp_),
fluid_(psf.fluid_),
cp_(psf.cp_,mapper),
rho_(psf.rho_,mapper),
lastTimeStep_(psf.lastTimeStep_)
{}

massTransferCoupledFvPatchScalarField::massTransferCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
):
mixedFvPatchScalarField(p, iF),
temperatureCoupledBase(patch(), dict),
pName_(dict.lookupOrDefault<word>("p", "p")),
UName_(dict.lookupOrDefault<word>("U", "U")),
rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
muName_(dict.lookupOrDefault<word>("mu", "thermo:mu")),
TnbrName_(dict.lookupOrDefault<word>("Tnbr", "T")),
specieName_(dict.lookupOrDefault<word>("specie", "none")),
sourceTerm_(dict.lookupOrDefault<word>("sourceTerm","eva")),
SolidDivCoe_(dict.lookupOrDefault<word>("solidCoeff","Ds")),
liquid_(nullptr),
liquidDict_(),
mass_(patch().size(), Zero),
massOld_(patch().size(), Zero),
Tvap_(0.0),
myKDelta_(patch().size(), Zero),
dmHfg_(patch().size(), Zero),
mpCpTp_(patch().size(), Zero),
Mcomp_(0.0),
fluid_(false),
cp_(patch().size(), Zero),
rho_(patch().size(), Zero),
lastTimeStep_(0.0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalIOErrorInFunction(dict)
            << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalIOError);
    }
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    if(fluid_)
    {
        dict.readEntry("carrierMolWeight",Mcomp_);
        dict.readEntry("Tvap",Tvap_);
        liquidDict_ = dict.subDict("liquid");
        liquid_ = liquidProperties::New(liquidDict_.subDict(specieName_));
    }
    lastTimeStep_ = patch().boundaryMesh().mesh().time().value();
}

massTransferCoupledFvPatchScalarField::massTransferCoupledFvPatchScalarField
(
    const massTransferCoupledFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
):
mixedFvPatchScalarField(psf, iF),
temperatureCoupledBase(patch(), psf),
TnbrName_(psf.TnbrName_),
specieName_(psf.specieName_),
sourceTerm_(psf.sourceTerm_),
SolidDivCoe_(psf.SolidDivCoe_),
liquid_(psf.liquid_.clone()),
liquidDict_(psf.liquidDict_),
mass_(psf.mass_),
massOld_(psf.massOld_),
myKDelta_(psf.myKDelta_),
dmHfg_(psf.dmHfg_),
mpCpTp_(psf.mpCpTp_),
Mcomp_(psf.Mcomp_),
fluid_(psf.fluid_),
lastTimeStep_(psf.lastTimeStep_)
{}

// * * * * * * * * * * * * * * * * * * * update * * * * * * * * * * * * * * * * * * * // 
void massTransferCoupledFvPatchScalarField::updateCoeffs()
{
    if(updated())
    {
        return;
    }

    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;
    //获取耦合边界信息
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch().patch());
    //计算单元的面积
    const scalarField& magSf = patch().magSf();
    //获取边界的网格
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    //获取区域的网格
    const polyMesh& nbrMesh = mpp.sampleMesh();
    //获取该边界网格的索引
    const label samplePatchi = mpp.samplePolyPatch().index();
    //获得对边界网格的引用
    const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];
    //获取温度场
    scalarField Tinternal(patchInternalField());

    scalarField& Tpatch = *this;
    //网格的内部场的引用
    const volScalarField& TmeshField = static_cast<const volScalarField&>
    (
        internalField()
    );
    //获取上一个时间步的场的数据
    scalarField TpatchOld
    (
        TmeshField.oldTime().boundaryField()[patch().index()]
    );

    //寻找要计算的场的数据，并存储到nbrTp中
    const fvPatchScalarField& nbrTp = nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);
    //获取固体上的物料质量分数
    //const fvPatchScalarField& YpatchS = nbrPatch.lookupPatchField<volScalarField, scalar>(specieName_);

    typedef massTransferCoupledFvPatchScalarField thisType;
    if (!isA<thisType>(nbrTp))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << patch().name() << " is of type " << thisType::typeName
            << endl << "The neighbouring patch field " << TnbrName_ << " on "
            << nbrPatch.name() << " is required to be the same, but is "
            << "currently of type " << nbrTp.type() << exit(FatalError);
    }
    const thisType& nbrField = refCast<const thisType>(nbrTp);

    //Swap to obtain full local values of neighbour internal field
    scalarField nbrInternalField(nbrField.patchInternalField());
    mpp.distribute(nbrInternalField);
    //温度的扩散率乘体心到面心的距离
    myKDelta_ = kappa(Tpatch)*patch().deltaCoeffs();
    //存储质量变化
    scalarField dm(patch().size(), Zero);
    //质量变化引起的能量的变化
    scalarField dHspec(patch().size(),Zero);
    //相变能
    scalarField hPhaseChange(patch().size(), Zero);
    //质量*相变能
    scalarField dmhPhaseChange(patch().size(), Zero);
    //质量所含能量
    scalarField hRemovedMass(patch().size(), Zero);
    //质量*质量所含能量
    scalarField dmhRemoveMass(patch().size(), Zero);

    if(fluid_)
    {
        const scalar dt = mesh.time().deltaTValue();
        //如果到了最终时刻，不再计算传质传热
        if(mesh.time().value() != lastTimeStep_)
        {
            lastTimeStep_ = mesh.time().value();
            massOld_ = mass_;
        }
        /*-----------------------------------------------------------------------------------
        //流体侧的通过thermo找到物性
        const rhoReactionThermo & thermo = mesh.lookupObject<rhoReactionThermo>(basicThermo::dictName);
        //获取组分的数据
        const basicSpecieMixture& composition = thermo.composition();
        //获取需要的组分的数据标签
        label specieIndex = composition.species()[specieName_];
        ----------------------------------------------------------------------------------*/
        //固体侧的定义上搜寻网格
        //const solidThermo& solidT = mesh.lookupObject<solidThermo>(basicThermo::dictName);
        //固体侧的密度
        //const volScalarField solidRho = solidT.rho();
        //固体侧边界的物料
        const fvPatchScalarField& YpatchS = nbrPatch.lookupPatchField<volScalarField,scalar>(specieName_);
        //固体侧边界的密度
        const fvPatchScalarField& rhoPatchS = nbrPatch.lookupPatchField<volScalarField,scalar>("rho");
        //固体侧边界的温度
        const fvPatchScalarField& TPatchS = nbrPatch.lookupPatchField<volScalarField,scalar>("T");
        //固体侧的温度场
        const volScalarField& Tsolid = nbrMesh.lookupObject<volScalarField>("T");
        //固体的扩散系数
        const fvPatchScalarField& DPatchS = nbrPatch.lookupPatchField<volScalarField,scalar>(SolidDivCoe_);
        //固体的边界旁的网格上的物料
        const scalarField YinternalS(YpatchS.patchInternalField());
        //固体侧的边界旁的网格上的密度
        const scalarField rhoInternalS(rhoPatchS.patchInternalField());
        //获取颗粒的网格体积
        const scalarField volumsS(nbrMesh.cellVolumes());
        //获取边界上的单元label
        const labelList patchLabelS(nbrPatch.faceCells());
        scalarField volInternalS(patch().size(),Zero);
        //提取边界单元的体积
        forAll(patchLabelS,cellI)
        {
            volInternalS[cellI] = volumsS[cellI];
        }
        //获取体心到面心的值
        const scalarField myDelta(patch().deltaCoeffs());
        //定义Cp
        scalarField cp(patch().size(), Zero);
        //定义密度
        scalarField liquidRho(patch().size(), Zero);


        //获取边界上的组分数据
        const fvPatchScalarField& Ypatch = patch().lookupPatchField<volScalarField, scalar>(specieName_);
        //边界的p
        const fvPatchScalarField& pPatch = patch().lookupPatchField<volScalarField, scalar>("p");
        //边界的密度
        const fvPatchScalarField& rhoPatch = patch().lookupPatchField<volScalarField, scalar>("rho");
        //动力粘度
        const fvPatchScalarField& muPatch = patch().lookupPatchField<volScalarField, scalar>("thermo:mu");
        //运动粘度
        const fvPatchScalarField& nutPatch = patch().lookupPatchField<volScalarField, scalar>("nut");        
        //提取边界上的数据
        //物料的
        scalarField Yinternal(Ypatch.patchInternalField());
        //p数据
        const scalarField pInternal(pPatch.patchInternalField());
        //密度数据
        const scalarField rhoInternal(rhoPatch.patchInternalField());
        //网格的体积
        const scalarField volumeInternel(mesh.cellVolumes());
        //获取单元的标签list
        const labelList& faceCells = patch().faceCells();

        forAll(Tpatch , faceI)
        {
            //面上的温度
            const scalar Tface = Tpatch[faceI];
            //网格温度
            const scalar Tcell = Tinternal[faceI];
            //面压强
            const scalar pFace = pPatch[faceI];
            //网格压强
            const scalar pCell = pInternal[faceI];
            //面动力粘度
            const scalar muFace = muPatch[faceI];
            //网格运动粘度
            const scalar nuCell = 2.74e-5;
            //面密度
            const scalar rhoFace = rhoPatch[faceI];
            //网格密度
            const scalar rhoCell = rhoInternal[faceI];
            //网格动力粘度
            const scalar mucell = nuCell * rhoCell;
            //湍流动力粘度
            // const scalar mutFace = nutPatch[faceI]*rhoFace;
            //面饱和蒸汽压
            const scalar pSatFace = liquid_->pv(pFace,Tcell);
            //网格饱和蒸汽压
            const scalar pSatCell = liquid_->pv(pCell,Tcell);
            //H2O的分子量
            const scalar Mv = liquid_->W();
            //物性的网格质量分数
            const scalar Ycell = Yinternal[faceI];
            //cell ---> face 的距离
            const scalar deltaFace = myDelta[faceI];
            //计算面的比热
            cp[faceI] = liquid_->Cp(pFace,Tface);
            
            /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\ 
            //计算面的汽化潜热
            hPhaseChange[faceI] = liquid_->hl(pFace,Tface);
            //计算面的蒸发质量
            hRemovedMass[faceI] = liquid_->Hs(specieIndex,pCell,Tcell);
            \* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
            //质量蒸发速率
            //dm[faceI] = volInternalS[faceI]*rhoPatchS[faceI]*YinternalS[faceI]*EvaporationRate(Tcell);
            //汽化热
            hPhaseChange[faceI] = liquid_->hl(pFace,Tface);
            //内能
            hRemovedMass[faceI] = liquid_->Hs(pFace,Tface);

            //能量传递的内能增加量 Q = cp*V*rho*Y*ΔT 
            scalar Q = cp[faceI]*volInternalS[faceI]*rhoPatchS[faceI]*YinternalS[faceI]*(Tsolid[faceI] - Tsolid.oldTime()[faceI]);
            //蒸发速率引起的内能变化 Q = k*m*H
            scalar Q_dot = EvaporationRate(Tcell) * volInternalS[faceI]*rhoPatchS[faceI]*YinternalS[faceI] * hPhaseChange[faceI];

            if (Q < Q_dot)
            {
                //如果能量传递的总值小于蒸发速率所需要的能量
                dm[faceI] = Q / hPhaseChange[faceI];
            }else{
                //能量传递的总值大于蒸发所需要的能量，按照蒸发速率来确定dm
                dm[faceI] = volInternalS[faceI]*rhoPatchS[faceI]*YinternalS[faceI]*EvaporationRate(Tcell);
            }

            //YpatchS[faceI] = YinternalS[faceI]/(1 - EvaporationRate(Tcell)*deltaFace/DPatchS[faceI]);
            //Yinternal[faceI] = dm[faceI]/2.2e-3/rhoPatch[faceI]/volumeInternel[faceI]*deltaFace - Ypatch[faceI];
        }
        //计算质量传输量
        mass_ = massOld_ + dm * dt ;
        mass_ = max(mass_, scalar(0));
        //计算蒸发损失的能量
        dmHfg_ = dm * hPhaseChange;
        //质量流出的能量
        dHspec = dm * hRemovedMass;


        //输出质量流量数据
        scalarField &massFluxOut = outputScalarField(
                specieName_+"MassFlux",
                dimMass/dimTime,
                refCast<const fvMesh>(mesh)
            ).boundaryFieldRef()[patch().index()];

        massFluxOut = dm;
        /*
        forAll(faceCells,faceI)
        {
            const label cellI = faceCells[faceI];

            //Ypatch[cellI] += dm[faceI] * ;
        }
        */
    }
    scalarField  mpCpTpNbr(patch().size(),Zero);
    scalarField  dmHfgNbr(patch().size(),Zero);
    //固体侧的处理
    if(!fluid_)
    {
        mpCpTpNbr = nbrField.mpCpTp_;
        mpp.distribute(mpCpTpNbr);

        dmHfgNbr = nbrField.dmHfg_;
        mpp.distribute(dmHfgNbr);        
    } 

    scalarField KDeltaNbr_;
    //计算温度变化引起的热量传递值
    KDeltaNbr_ = nbrField.kappa(nbrField) * nbrPatch.deltaCoeffs();
    mpp.distribute(KDeltaNbr_);

    const scalarField dmHfg(dmHfgNbr + dmHfg_);
    const scalarField mpCpTp(mpCpTp_ + mpCpTp_);

    scalarField alpha(KDeltaNbr_ + mpCpTp);

    valueFraction() = alpha/(alpha + myKDelta_);
    refValue() = (KDeltaNbr_*nbrInternalField + mpCpTp*TpatchOld + dmHfg) / alpha ;

    if (fluid_)
    {
        scalar Qdm = gSum(dm);
        scalar QMass = gSum(mass_);
        scalar Qt = gSum(myKDelta_ * (Tpatch - Tinternal) * magSf);
        scalar QtSolid = gSum(KDeltaNbr_*(Tpatch - nbrInternalField)*magSf);

        scalar Q = gSum(kappa(Tpatch)*patch().magSf()*snGrad());
        scalar Qhgf = gSum(dmHfg_*magSf);
        scalar Qspec = gSum(dHspec*magSf);

        Info<< mesh.name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :" << nl
            << " heat transfer rate:" << Q << endl
            << "    Total mass flux [Kg/s]:                            "
            << Qdm << nl
            << "    Total mass on the wall [Kg]:                       "
            << QMass << nl
            << "    Total heat (>0 leaving the wall to the fluid) [W]: "
            << Qt << nl
            << "    Total heat (>0 leaving the wall to the solid) [W]: "
            << QtSolid << nl
            << "    Total latent heat released [W]:                    "
            << Qhgf << nl
            << "    Total specific heat removed from fluid [W]:        "
            << Qspec << nl
            << " walltemperature "
            << " min:" << gMin(Tpatch)
            << " max:" << gMax(Tpatch)
            << " avg:" << gAverage(Tpatch)
            << endl;

    }

    UPstream::msgType() = oldTag ;

}


void massTransferCoupledFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;

    if (fluid_)
    {
        os.writeKeyword("specie") << specieName_ << token::END_STATEMENT << nl;
        os.writeKeyword("carrierMolWeight") << Mcomp_ << token::END_STATEMENT << nl;
        mass_.writeEntry("mass", os);
        word liq = "liquid";
        os << token::TAB << token::TAB << liq;
        liquidDict_.write(os);
    }

    temperatureCoupledBase::write(os);
}

makePatchTypeField
(
    fvPatchScalarField,
    massTransferCoupledFvPatchScalarField
);


}//end namespace compressible
}//end namespace foam