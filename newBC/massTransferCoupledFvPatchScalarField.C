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
{/*
    if(T < scalar(430))
    {
        scalar A = 5.13e10;//Amoi   /s
        scalar E = 27360;//Emoi     J/mol
        scalar R = 8.31446261815;//R     J/kmol
        scalar k = exp(-1 * E / R / T);
        return k; 
    }else{
        return scalar(0.9);
    }*/
    scalar k = 0.24261 * exp( -2530.2 / T );
    return k;
}
//计算舍伍德数
scalar massTransferCoupledFvPatchScalarField::Shnumber(scalar Re,scalar Sc) const
{
    scalar shnumber_ = 2.0+0.552*sqrt(Re)*cbrt(Sc);
    return shnumber_;
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

/* * * * * * * * * * * * 构造函数 * * * * * * * * * * * */

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
sourceTerm_("HSource"),
L_(0.001),
SolidDivCoe_("Dsolid"),
liquid_(nullptr),
liquidDict_(nullptr),
mass_(patch().size(), Zero),
massOld_(patch().size(), Zero),
myKDelta_(patch().size(), Zero),
dmHfg_(patch().size(), Zero),
mpCpTp_(patch().size(), Zero),
Mcomp_(0.0),
fluid_(false),
cp_(patch().size(), Zero),
rho_(patch().size(), Zero),
thickness_(patch().size(), Zero)
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
L_(psf.L_),
SolidDivCoe_(psf.SolidDivCoe_),
liquid_(psf.liquid_.clone()),
liquidDict_(psf.liquidDict_),
mass_(psf.mass_,mapper),
massOld_(psf.massOld_,mapper),
myKDelta_(psf.myKDelta_,mapper),
dmHfg_(psf.dmHfg_,mapper),
mpCpTp_(psf.mpCpTp_,mapper),
Mcomp_(psf.Mcomp_),
fluid_(psf.fluid_),
cp_(psf.cp_,mapper),
rho_(psf.rho_,mapper),
thickness_(psf.thickness_, mapper)
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
specieName_(dict.lookupOrDefault<word>("specie", "H2O")),
sourceTerm_(dict.lookupOrDefault<word>("sourceTerm","HSource")),
L_(0.001),
SolidDivCoe_(dict.lookupOrDefault<word>("solidCoeff","Dsolid")),
liquid_(nullptr),
liquidDict_(),
mass_(patch().size(), Zero),
massOld_(patch().size(), Zero),
myKDelta_(patch().size(), Zero),
dmHfg_(patch().size(), Zero),
mpCpTp_(patch().size(), Zero),
Mcomp_(0.0),
fluid_(false),
cp_(patch().size(), Zero),
rho_(patch().size(), Zero),
thickness_(patch().size(), Zero)
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

    if (dict.found("refValue"))
    {
        // Full restart
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
        fluid_ = true;
        Info << "************************ found fluid side ************************" << endl;
        //获取边界的网格信息
        // const polyMesh& mesh = patch().boundaryMesh().mesh();
        // //查找当前网格的网格名；
        // const fvPatch& yp = refCast<const fvMesh>(mesh).boundary()[patch().index()];
        // const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch().patch());
        // const fvPatch& sp = refCast<const fvMesh>(mpp.sampleMesh()).boundary()[mpp.samplePolyPatch().index()];
        // Info <<"\tmapped patch name is : " << yp.name() << "\n  and it's nbr patch name is : "<<sp.name()<< endl;

        /*forAll(yp,faceI){
            Info << faceI << "\t" ;
        }Info << endl;
        forAll(sp ,faceI){
            Info << faceI << "\t" ;
        }Info << endl;*/
    }else
    {
        fluid_ = false;
        Info << "************************ found solid side ************************" << endl;
        //获取边界的网格信息
        // const polyMesh& mesh = patch().boundaryMesh().mesh();
        // //查找当前网格的网格名；
        // const fvPatch& yp = refCast<const fvMesh>(mesh).boundary()[patch().index()];
        // const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch().patch());
        // const fvPatch& sp = refCast<const fvMesh>(mpp.sampleMesh()).boundary()[mpp.samplePolyPatch().index()];
        // Info <<"\tmapped patch name is : " << yp.name() << "\n  and it's nbr patch name is : "<<sp.name()<< endl;
    }

    if(fluid_)
    {
        dict.readEntry("carrierMolWeight",Mcomp_);
        dict.readEntry("CharaLength",L_);
        liquidDict_ = dict.subDict("liquid");
        liquid_ = liquidProperties::New(liquidDict_.subDict(specieName_));

        if (dict.found("thickness"))
        {
            scalarField& Tp = *this;
            const scalarField& magSf = patch().magSf();
            // Assume initially standard pressure for rho calculation
            scalar pf = 1e5;
            thickness_ = scalarField("thickness", dict, p.size());
            forAll(thickness_, i)
            {
                mass_[i] = thickness_[i]*liquid_->rho(pf, Tp[i])*magSf[i];
                massOld_[i] = mass_[i];
            }
        }
    }
}

massTransferCoupledFvPatchScalarField::massTransferCoupledFvPatchScalarField
(
    const massTransferCoupledFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
):
mixedFvPatchScalarField(psf, iF),
temperatureCoupledBase(patch(), psf),
pName_(psf.pName_),
UName_(psf.UName_),
rhoName_(psf.rhoName_),
muName_(psf.muName_),
TnbrName_(psf.TnbrName_),
specieName_(psf.specieName_),
sourceTerm_(psf.sourceTerm_),
L_(psf.L_),
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
cp_(psf.cp_),
rho_(psf.rho_),
thickness_(psf.thickness_)
{}

// * * * * * * * * * * * * * * * *  mapper function * * * * * * * * * * * * * * * * * // 

void massTransferCoupledFvPatchScalarField::autoMap(const fvPatchFieldMapper& m)
{
     mixedFvPatchScalarField::autoMap(m);

     if (fluid_)
     {
         mass_.autoMap(m);
         myKDelta_.autoMap(m);
         dmHfg_.autoMap(m);
         mpCpTp_.autoMap(m);
         rho_.autoMap(m);
     }
 }

void massTransferCoupledFvPatchScalarField::rmap(const fvPatchScalarField& ptf,const labelList& addr)
{
     mixedFvPatchScalarField::rmap(ptf,addr);

     const massTransferCoupledFvPatchScalarField& tiptf = refCast<const massTransferCoupledFvPatchScalarField>(ptf);

     if (fluid_)
     {
         mass_.rmap(tiptf.mass_,addr);
         myKDelta_.rmap(tiptf.myKDelta_,addr);
         dmHfg_.rmap(tiptf.dmHfg_,addr);
         mpCpTp_.rmap(tiptf.mpCpTp_,addr);
         rho_.rmap(tiptf.rho_,addr);
     }
 }


// * * * * * * * * * * * * * * * * * * * update * * * * * * * * * * * * * * * * * * * // 
void massTransferCoupledFvPatchScalarField::updateCoeffs()
{
    if(updated())
    {
        return;
    }
	
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;
    
	if (fluid_)
    {
        Pout << "this is fluid regions properties: \n" ;
    }else{
        Pout << "this is solid regions properties: \n" ;
    }
    //换名字
    typedef massTransferCoupledFvPatchScalarField thisType;
    //获取耦合边界信息
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(patch().patch());
    //获取边界的网格信息
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    //获取邻区的边界的网格信息
    const polyMesh& nbrMesh = mpp.sampleMesh();
    //获取当前伶区边界的网格index
    const label samplePatchi = mpp.samplePolyPatch().index();
    //获取当前邻近区域的边界的patch信息
    const fvPatch& nbrPatch = refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];
	//获取邻区的温度场
    const thisType& nbrField = refCast<const thisType>(nbrPatch.lookupPatchField<volScalarField,scalar>(TnbrName_));
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);
    //获取当前区域的边界的网格体心值
    scalarField Tinternal(patchInternalField());
    //获取当前区域的面心上的数据
    scalarField& Tp = *this;
    //查找当前边界的边界信息提取
    const fvPatch& yp = refCast<const fvMesh>(mesh).boundary()[patch().index()];
    //获取所有边界的单元面积
    const scalarField& magSf = yp.magSf();
    //获取内部的场的信息到T中
    const volScalarField& T = static_cast<const volScalarField&>(internalField());
    //获取边界上场的数值
    scalarField Tf(T.boundaryField()[patch().index()]);
    //获取上一个该边界的时间步的信息
    scalarField TOld(T.oldTime().boundaryField()[patch().index()]);
    //体心数据
    scalarField Tin(patchInternalField());
    //计算此边界的kappa
    scalarField K(this->kappa(*this));
    //获取邻区的htc
    scalarField nbrK(nbrField.kappa(*this));
    //Pout << "line 1 " << nl;
	if (fluid_)
    {
        //Pout << "fluid side calculate the tempNbrK " << nl;
		const fvPatchScalarField& tempnbrk = nbrPatch.lookupPatchField<volScalarField,scalar>("solidKappa");
        nbrK = tempnbrk.internalField();
    }else if(yp.size() > 0)
    {
		//Pout << "solid side calculate the tempNbrK " << nl;
        const fvPatchScalarField& tempnbrk = yp.lookupPatchField<volScalarField,scalar>("solidKappa");
        //Pout << "the patch's temp nbrk size is :" << tempnbrk.size() << nl;
		K = tempnbrk.internalField();
		//Pout << "and the patch's size is : " << patch().size() << nl;
    }
    mpp.distribute(nbrK);
    //获取deltaCoeffs，体心到面心距离
    //Pout << "nbrDeltaCoeffs update " << nl;
	scalarField nbrDeltaCoeffs(nbrPatch.deltaCoeffs());
    mpp.distribute(nbrDeltaCoeffs);
    //计算 k*Δ htc
	//Pout << "KDeltaNbr coeffs update " << nl;
    scalarField KDeltaNbr(nbrK*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);
    //计算此边界的Kdelta
	//Pout << "calculate this patch KDelta " << nl;
    myKDelta_ = K * patch().deltaCoeffs();
	//Pout << "line 2 " << nl;
    //检测输出信息
    /*
    Info << "============================================================================================================\n"
         << "this patch name:"<<yp.name() <<", size is :"<< yp.size()<<", T[200]="<<T[200]<<", Told[200]="<<TOld[200]
         << ", Tf's size:" << Tf.size() << ", TOld's size:" <<TOld.size()
         << "\nthis mapped patch name is " << nbrPatch.name() << ", size is :" << nbrPatch.size() << " and nbrK size:" << nbrK.size()
         << ", nbrK[200]:" << nbrK[200] <<", myKDelta_ size:" <<myKDelta_.size() << ", myKDelta_[200]:" << myKDelta_[200]
         << ", nbrDeltaCoeffs size:"<< nbrDeltaCoeffs.size() <<", nbrDeltaCoeffs[200]:"<<nbrDeltaCoeffs[200]
         << ", KDeltaNbr size:"<<KDeltaNbr.size()<<", KDeltaNbr[200]:" << KDeltaNbr[200] 
         << "\n============================================================================================================"<<endl;
    */
    /*-------------------------------------------------------------------------------------*\
                                        计算存储的数据值
    \*-------------------------------------------------------------------------------------*/
    //存储的质量流量kg/s/m2
    scalarField dm(yp.size(),Zero);

    if(fluid_)
    {
        //获取时间步长
        const scalar dt = mesh.time().deltaTValue();
        /*-------------------------------------------------------------------*\
                                固体侧 solid
        \*-------------------------------------------------------------------*/
        //获取物料边界
        fixedGradientFvPatchField<scalar>& YpatchSolid = const_cast<fixedGradientFvPatchField<scalar>&>
        (refCast<const fixedGradientFvPatchField<scalar>>(nbrPatch.lookupPatchField<volScalarField,scalar>(specieName_)));
        //获取物料体心数据
        const scalarField Yisolid(YpatchSolid.patchInternalField());
        //获取边界的网格体积
        //const fvPatchScalarField& VolSolidPatchBoundary = nbrPatch.lookupPatchField<volScalarField,scalar>("cellVolume");
        //const scalarField& VolSolidPatch = VolSolidPatchBoundary.patchInternalField();
        //获取边界的源项
        fvPatchScalarField& HsourceSolidpatch = const_cast<fvPatchScalarField&>
        (refCast<const fvPatchScalarField>(nbrPatch.lookupPatchField<volScalarField,scalar>(sourceTerm_)));
        //获取固体的密度
        const fvPatchScalarField& rhoSolidpatch = nbrPatch.lookupPatchField<volScalarField,scalar>("thermo:rho");
        //获取固体侧的扩散率
        const fvPatchScalarField& Dsolid = nbrPatch.lookupPatchField<volScalarField,scalar>(SolidDivCoe_);
        scalarField Dsolidinter(Dsolid.internalField());
        //获取固体上的温度值和旧值
        const volScalarField& TSolid = static_cast<const volScalarField&>(nbrField.internalField());
        scalarField TpatchSolidOld(TSolid.internalField());
        //patch上的面积
        scalarField magSfSolid(nbrPatch.magSf());
        scalarField soliddelta(nbrPatch.deltaCoeffs());
        /*
        Info << "solid properties : solidpatch cellVolume[100]:" << magSfSolid[200]/soliddelta[200]
             << ", magSfSolid[100]:"<< magSfSolid[200]
             << ", soliddelta[100]:"<< soliddelta[200]
             << ", HsourceSolid[200]:"<< HsourceSolidpatch[200]
             << ", Dsolid[200]:" << Dsolid[200]
             << ", solidTemp[200]:" << nbrField[200]
             << ", TpatchSolidOld[200]:" << TpatchSolidOld[200]
             << " \n" << endl;
        */
        /*-----------------------------------------------------------------*\
                                气体侧 gas
        \*-----------------------------------------------------------------*/
        //获取物料边界
        fixedGradientFvPatchField<scalar>& YpatchFluid = const_cast<fixedGradientFvPatchField<scalar>&>
        (refCast<const fixedGradientFvPatchField<scalar>>(yp.lookupPatchField<volScalarField,scalar>(specieName_))); 
        //获取物料的体心数据
        const scalarField Yifluid(YpatchFluid.patchInternalField());
        //获得流体域的压力场
        const fvPatchScalarField& ppatchFluid = yp.lookupPatchField<volScalarField,scalar>(pName_);
        //获得流体域的边界温度
        const fvPatchScalarField& TpatchFluid = yp.lookupPatchField<volScalarField,scalar>(TnbrName_);
        //获得流体的流速
        const fvPatchVectorField& UpatchFluid = yp.lookupPatchField<volVectorField,vector>(UName_);
        const vectorField UFluid(UpatchFluid.patchInternalField()); 
        //获得边界体心的温度
        scalarField TinPatchFluid(TpatchFluid.patchInternalField());
        //获得边界的密度
        const fvPatchScalarField& rhoPatchfluid = yp.lookupPatchField<volScalarField,scalar>(rhoName_);
        //获取边界的运动粘度
        const fvPatchScalarField& muPatchfluid = yp.lookupPatchField<volScalarField,scalar>(muName_);

        //存储比热
        scalarField cp(yp.size(),Zero);
        //存储潜热
        scalarField hfg(yp.size(),Zero);
        //流体密度
        scalarField liquidRho(yp.size(),Zero);
        /*
        Info << "fluid properties : fluid TpatchFluid[100]:" << TpatchFluid[100]
             << ", ppatchFluid[200]:"<< ppatchFluid[200]
             << " \n" << endl;   */
        /*-----------------------------------------------------------------*\
                                计算 calculate
        \*-----------------------------------------------------------------*/
        scalarField yvp(yp.size(),Zero);
		scalarField yvp_s(yp.size(),Zero);
        //Pout << " ======= Start calculate of patch =======" << endl;     
        forAll(yp,faceI)
        {
            // Info << "nbrK : " << nbrK[faceI] << ", K : " << K[faceI] << ", Yisolidtmp : " << Yisolid[faceI] <<" , YpatchSolid:"<< YpatchSolid[faceI] << endl;
            // scalar TFfluid = Tp[faceI];
            // Info << "TFfluid" <<TFfluid<< endl;
            scalar Tinfluid = TinPatchFluid[faceI];
            // Info << "Tinfluid" << Tinfluid<< endl;
            scalar Tinsolid = nbrField[faceI];
            // Info << "Tinsolid" << Tinsolid<< endl;
            scalar pFfluid = ppatchFluid[faceI];
            // Info << "pFfluid" << pFfluid<< endl;
            scalar rhofluid = rhoPatchfluid[faceI];
            // Info << "rhofluid" << rhofluid<< endl;
            scalar rhosolid = rhoSolidpatch[faceI];
            // Info << "rhosolid" << rhosolid<< endl;
            scalar Yi = Yifluid[faceI];
            // Info << "Yi" << Yi<< endl;
            scalar Dsolidtmp = Dsolidinter[faceI];
            // Info << "Dsolidtmp" << Dsolidtmp<< endl;
            scalar Yisolidtmp = Yisolid[faceI];
            // Info << "Yisolidtmp" << Yisolidtmp<< endl;
            liquidRho[faceI] = rhofluid;
            // Info << "liquidRho[faceI]" << liquidRho[faceI]<< endl;
            cp[faceI] = liquid_->Cp(pFfluid,Tinfluid);
            // Info << "cp[faceI]" << cp[faceI]<< endl;
            hfg[faceI] = liquid_->hl(pFfluid,Tinfluid);
            // Info << "hfg[faceI]" << hfg[faceI]<< endl;
            scalar DeltaSolid = nbrDeltaCoeffs[faceI];
            //动力粘度计算
            scalar nuf = muPatchfluid[faceI]/rhofluid;            
            // Info << "nuf, pFfluid: "<< pFfluid <<", TFfluid: "<< TFfluid << endl;
            //蒸汽扩散系数计算
            const scalar Dpatch = liquid_->D(pFfluid,Tinfluid);            
            // Info << "Dpatch" << endl;
            //获得蒸汽的分子质量
            const scalar Mv = liquid_->W();            
            // Info << "Mv" << endl;
            //施密特数
            scalar Sc = nuf / Dpatch / rhofluid ;            
            // Info << "Sc" << endl;
            const vector Ufluid = UFluid[faceI];            
            // Info << "Ufluid" << endl;
            //雷诺数计算
            scalar Re = mag(Ufluid) * L_ / nuf;            
            // Info << "Re" << endl;
            //蒸汽的当前状态的饱和蒸汽压
            const scalar Psat = liquid_->pv(pFfluid,Tinfluid);            
            // Info << "Psat" << endl;
            //RH计算
            const scalar invMwmean = Yi/Mv + (1 - Yi) / Mcomp_;            
            // Info << "invMwmean" << endl;
            scalar Xv = Yi/Mv/invMwmean;            
            // Info << "Xv" << endl;
            scalar RH = min(pFfluid * Yi /(Yi + 0.622) / Psat,0.998);
            // Info << "RH" << endl;
            scalar Ysat = 0.622 * Psat / (pFfluid - Psat);
            scalar Weq = 0.295-0.045*Foam::log(-(Tinsolid-238.14)*Foam::log(RH));

            scalar yvp_tmp = (Yisolidtmp - Weq)*DeltaSolid;
            scalar dmtemp = yvp_tmp * Dsolidtmp * rhosolid;

            // Info << "5 RH:" << RH << ", Re: " << Re << ", Dpatch: " << Dpatch << ", Dsolid: " << Dsolidtmp<<endl;
            //计算当前的传质系数 m/s
            const scalar Hm = Dpatch * Shnumber(Re,Sc) / L_;
            //单位面积的质量流量 kg/m^2/s
            scalar dmfluid = rhofluid * Hm * max((Ysat - Yi),0.0);
            scalar dmsolid = rhosolid * Hm * (Yisolidtmp - Weq);
            dm[faceI] = min(min(dmfluid,dmsolid),dmtemp);
            //传质量 kg
            mass_[faceI] = dm[faceI] * magSf[faceI] * dt;
            // Info << "7 Hm: "<<Hm << endl;
            //蒸发量的能量      J
            const scalar q = dm[faceI] * dt * magSf[faceI] * hfg[faceI];
            //梯度计算
            yvp[faceI] = dm[faceI]/ Dpatch/ rhofluid;
			yvp_s[faceI] = scalar(-1.0) * dm[faceI] / Dsolidtmp / rhosolid;
        }
        //Pout << " ======= End calculate of patch =======" << endl;
        //传递梯度值，流体为正，固体为负，产生传递现象
        YpatchFluid.gradient() = yvp;
        YpatchSolid.gradient() = yvp_s;
        //更正总的传递质量
        mass_ = max(mass_,scalar(0));
        //输出当前累加的传递质量
        scalarField& massTotadd = outputScalarField(
                specieName_ + "massTot",
                dimMass,
                refCast<const fvMesh>(mesh)
            ).boundaryFieldRef()[patch().index()];
        massTotadd = mass_;
        mpCpTp_ = dm * cp;
        dmHfg_ = dm * hfg;
    }
    
    scalarField  mpCpTpNbr(yp.size(),Zero);
    scalarField  dmHfgNbr(yp.size(),Zero);
    //固体侧的处理
	//Pout << "calculate the solid side ! " << nl;
    if(!fluid_)
    {
        mpCpTpNbr = nbrField.mpCpTp_;
        mpp.distribute(mpCpTpNbr);

        dmHfgNbr = nbrField.dmHfg_;
        mpp.distribute(dmHfgNbr);        
    } 

    const scalarField dmHfg(dmHfgNbr - dmHfg_);
    const scalarField mpCpdt(mpCpTpNbr - mpCpTp_);

    scalarField alpha(KDeltaNbr - mpCpdt);

    valueFraction() = alpha/(KDeltaNbr + myKDelta_);

    refValue() = (KDeltaNbr * nbrIntFld - dmHfg - mpCpdt * TOld)/alpha;

    mixedFvPatchScalarField::updateCoeffs();
    
    if (fluid_ && debug)
    {
        scalar Qdm = gSum(dm);
        scalar QMass = gSum(mass_);
        scalar QtFluid = gSum(myKDelta_ * (Tp - Tin) * magSf);
        scalar QtSolid = gSum(KDeltaNbr * (nbrField - nbrIntFld) * magSf);
        scalar Q = gSum(kappa(Tp) * patch().magSf()*snGrad());
        scalar Qhgf = gSum(dmHfg_*magSf);
        Info << mesh.name() << ":" << patch().name() << ":" 
             << internalField().name() << "->"
             << nbrMesh.name() << ":" << nbrPatch.name() << ":"
             << internalField().name() << ":" << nl
             << " \t Total mass flux [kg/s] :" << Qdm << nl
             << " \t Total mass on the wall [Kg] :" << QMass << nl
             << " \t Total heat (>0 leaving the wall to the fluid) [W] :" << QtFluid << nl
             << " \t Total heat (>0 leaving the wall to the solid) [W] :" << QtSolid << nl
             << " \t the heat transfer rate is :" << Q << nl
             << " \t the hfg energy is :" << Qhgf << nl
             << " \t wall temperature :" 
             << " \t min:" << gMin(*this)
             << " \t max:" << gMax(*this)
             << " \t avg:" << gAverage(*this)
             << " \n " << endl;
    }

    UPstream::msgType() = oldTag ;
}


void massTransferCoupledFvPatchScalarField::write(Ostream& os) const
{
    mixedFvPatchScalarField::write(os);
    // os.writeEntryIfDifferent<word>("p", "p", pName_);
    // os.writeEntryIfDifferent<word>("U", "U", UName_);
    // os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    // os.writeEntryIfDifferent<word>("mu", "thermo:mu", muName_);
    // os.writeEntryIfDifferent<word>("Tnbr", "T", TnbrName_);

    if (fluid_)
    {
        os.writeEntryIfDifferent<word>("specie", "none", specieName_);

        os.writeEntry("carrierMolWeight", Mcomp_);

        os.writeEntry("CharaLength", L_);
        os.writeEntry("fluid", fluid_);
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
