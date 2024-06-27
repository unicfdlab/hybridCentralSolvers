/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
       hybridCentralSolvers | Copyright (C) 2016-2021 ISP RAS (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.
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

#include "vofTwoPhaseCentralFoam.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "totalPressureFvPatchScalarField.H"

namespace Foam
{
    defineTypeNameAndDebug(vofTwoPhaseCentralFoam, 0);
}

Foam::vofTwoPhaseCentralFoam::vofTwoPhaseCentralFoam(const fvMesh& mesh, pimpleControl& ctrl)
:

    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    pimple_ (ctrl),

    U_
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    phi_
    (
        "phi",
        (fvc::interpolate(U_))&mesh.Sf()
    ),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    p_rgh_
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_,
        p_rghPatchTypes()
    ),

    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    volumeFraction1_
    (
        IOobject
        (
            "volumeFraction1",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    volumeFraction2_
    (
        "volumeFraction2",
        1.0 - volumeFraction1_
    ),

    mixture_model_(*this, volumeFraction1_, volumeFraction2_, p_, T_),

    turbulence_
    (
        incompressible::turbulenceModel::New
         (
             U_,
             phi_,
             mixture_model_
         )
    ),

    dotVF1_
    (
        "dotVF1", 0.0*volumeFraction1_/mesh.time().deltaT()
    ),

    dotVF2_
    (
        "dotVF2", 0.0*volumeFraction2_/mesh.time().deltaT()
    ),

    phiVF1_
    (
        "phiVF1", 0.0*phi_
    ),

    phiVF2_
    (
        "phiVF2", 0.0*phi_
    ),

    interface_
    (
        volumeFraction1_,
        U_,
        *this
    ),

    C_
    (
        "C",
        sqrt(mixture_model_.gamma1()*mixture_model_.R1()*T_)
    ),

    Lambda_
    (
        "Lambda",
        0*volumeFraction2_
    ),

    HbyA_
    (
        IOobject
        (
            "HbyA",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 1, -1, 0, 0, 0, 0)
    ),

    oneByA_
    (
        IOobject
        (
            "rbyA",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(-1, 3, 1, 0, 0, 0, 0)
    ),

    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),

    onemkappa_(1.0 - kappa_),

// /***********************Tadmor-Kurganov Scheme*******************************/

    v_zero_
    (
        "v_zero",
        dimVolume/dimTime,
        0.0
    ),

    own_
    (
        IOobject
        (
            "own",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("own", dimless, 1.0)
    ),

    nei_
    (
        IOobject
        (
            "nei",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("nei", dimless, -1.0)
    ),

    phi1_own_
    (
        "phi1_own",
        phi_*mixture_model_.rho01()*0.0
    ),

    phi1_nei_
    (
        "phi1_nei",
        phi_*mixture_model_.rho01()*0.0
    ),

    phi2_own_
    (
        "phi2_own",
        phi_*mixture_model_.rho02()*0.0
    ),

    phi2_nei_
    (
        "phi2_nei",
        phi_*mixture_model_.rho02()*0.0
    ),

    rho1_own_
    (
        "rho1_own",
        fvc::interpolate(mixture_model_.rho1(), own_, "reconstruct(rho)")
    ),

    rho1_nei_
    (
        "rho1_nei",
        fvc::interpolate(mixture_model_.rho1(), nei_, "reconstruct(rho)")
    ),

    rho2_own_
    (
        "rho2_own",
        fvc::interpolate(mixture_model_.rho2(), own_, "reconstruct(rho)")
    ),

    rho2_nei_
    (
        "rho2_nei",
        fvc::interpolate(mixture_model_.rho2(), nei_, "reconstruct(rho)")
    ),

    alpha_own_
    (
        "alpha_own",
        own_
    ),

    alpha_nei_
    (
        "alpha_nei",
        1.0 - alpha_own_
    ),

    aSf_
    (
        "aSf",
        (fvc::interpolate(U_, own_, "reconstruct(U)")) & mesh.Sf()
    ),

    Cf_own_
    (
        "Cf_own",
        0.0*aSf_/mesh.magSf()
    ),

    Cf_nei_
    (
        "Cf_nei",
        Cf_own_
    ),

    CfSf_own_
    (
        "CfSf_own",
        0.0*aSf_
    ),

    CfSf_nei_
    (
        "CfSf_nei",
        CfSf_own_
    ),

    amaxSf_
    (
        "amaxSf",
        0.0*aSf_
    ),

    aphiv_own_
    (
        "aphiv_own",
        0.0*aSf_
    ),

    aphiv_nei_
    (
        "aphiv_nei",
        aphiv_own_
    ),

    phiv_own_
    (
        "phiv_own",
        0.0*aSf_
    ),

    phiv_nei_
    (
        "phiv_nei",
        phiv_own_
    ),

    phiHbyA_
    (
        "phiHbyA",
        phi_*0.0
    ),

    phiHbyA_own_
    (
        "phiHbyA_own",
        0.0*aSf_
    ),

    phiHbyA_nei_
    (
        "phiHbyA_nei",
        phiHbyA_own_
    ),

    rAUf_own_
    (
        "rAUf_own",
        0.0*linearInterpolate(oneByA_)
    ),
    rAUf_nei_
    (
        "rAUf_nei",
        rAUf_own_
    ),

    TSource1_
    (
        fvc::ddt(p_)
    ),

    TSource2_
    (
        fvc::ddt(p_)
    ),

    E1_
    (
        fvc::ddt(mixture_model_.rho1())
    ),

    E2_
    (
        fvc::ddt(mixture_model_.rho2())
    ),

    E_
    (
        0*E1_
    ),

    vF1face_
    (
        fvc::interpolate(volumeFraction1_, "reconstruct(volumeFraction1)")
    ),

    vF2face_
    (
        1.0 - vF1face_
    ),

    Q_
    (
        0.5*magSqr(U_)
    ),

    dpdt_
    (
        0.0*p_/p_.mesh().time().deltaT()
    ),
    
    gh_
    (
        IOobject
        (
            "gh",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0, 0, 0)
    ),

    ghf_
    (
        IOobject
        (
            "ghf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 2, -2, 0, 0, 0, 0)
    ),

    Wp_
    (
        "Wp_",
        1.0/(1.0 - (volumeFraction1_*mixture_model_.psi1() + volumeFraction2_*mixture_model_.psi2())*gh_)
    ),

    B_
    (
        "B",
        0.0*mixture_model_.rho01()*U_/U_.mesh().time().deltaT()
    ),

    phib_
    (
        "phib",
        0.0*phi_
    )

{
    HbyA_.primitiveFieldRef() = vector::zero;
    HbyA_.boundaryFieldRef() = vector::zero;
    oneByA_.primitiveFieldRef() = 1.0;
    oneByA_.boundaryFieldRef() = 1.0;
    Info<< "\nAll fields were created\n" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vofTwoPhaseCentralFoam::~vofTwoPhaseCentralFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vofTwoPhaseCentralFoam::saveOld()
{
    volumeFraction1_.oldTime();
    volumeFraction2_.oldTime();
    U_.oldTime();
    T_.oldTime();
    p_.oldTime();
    Q_.oldTime();
    p_rgh_.oldTime();
    phi_.oldTime();
    mixture_model_.saveOldTime();
}

void Foam::vofTwoPhaseCentralFoam::CharacteristicCourant()
{
    const auto &magSf = phi_.mesh().magSf();
    const auto &deltaCoeffs = phi_.mesh().deltaCoeffs();
    const auto &deltaT = phi_.mesh().time().deltaT();

    surfaceScalarField CfSf (max(CfSf_own_, CfSf_nei_));
    amaxSf_ = max(mag(phi_ + CfSf), mag(phi_ - CfSf));
    // surfaceScalarField uPlusC_pos =
    //     max(max(phi_own_ + CfSf_own_,phi_nei_ + CfSf_nei_),v_zero_)/magSf;
    // surfaceScalarField uPlusC_neg =
    //     min(min(phi_own_ - CfSf_own_,phi_nei_ - CfSf_nei_),v_zero_)/magSf;
    // surfaceScalarField uPlusC_max = max(uPlusC_pos,-uPlusC_neg);
    
    surfaceScalarField CCof (amaxSf_ * deltaCoeffs * deltaT / magSf);
    Info<< "max/min CCof:"
        << gMax(CCof)
        << "/"
        << gMin(CCof)
        << ", CCof dims = " << CCof.dimensions()
        << endl;
}

Foam::scalar Foam::vofTwoPhaseCentralFoam::FlowCourant()
{
    const auto& V = phi_.mesh().V().field();
    const auto& deltaT = phi_.mesh().time().deltaTValue();

    scalarField sumPhi
    (
        fvc::surfaceSum(mag(phi_))().primitiveField()
    );

    scalar CoMax = gMax((sumPhi)/V)*deltaT;
    return CoMax;
}

void Foam::vofTwoPhaseCentralFoam::massError1()
{
    auto& rho1 = mixture_model_.rho1();
    E1_ =
    (
        fvc::ddt(rho1) + fvc::div(phi1_own_ + phi1_nei_)
    );

    Info<< " max: rho1 " << max(rho1).value()
        << " min: " << min(rho1).value()
        << endl;
}


void Foam::vofTwoPhaseCentralFoam::massError2()
{
    auto& rho2 = mixture_model_.rho2();
    E2_ =
    (
        fvc::ddt(rho2) + fvc::div(phi2_own_ + phi2_nei_)
    );

    Info<< " max: rho2 " << max(rho2).value()
        << " min: " << min(rho2).value()
        << endl;
}


void Foam::vofTwoPhaseCentralFoam::TSource()
{
    dpdt_ = (p_ - p_.oldTime()) / p_.mesh().time().deltaT();
    Q_ = 0.5*magSqr(U_);
    Q_.rename("Q");

    const auto& rho1 = mixture_model_.rho1();
    const auto& rho2 = mixture_model_.rho2();
    TSource1_ =
    (
        fvc::ddt(rho1,Q_)
      + fvc::div(phi1_own_,Q_) + fvc::div(phi1_nei_,Q_)
      - dpdt_
      - fvc::Sp(E1_,Q_)
    );

    TSource2_ =
    (
        fvc::ddt(rho2,Q_)
      + fvc::div(phi2_own_,Q_) + fvc::div(phi2_nei_,Q_)
      - dpdt_
      - fvc::Sp(E2_,Q_)
    );
}

//* * * * * * * * * * * * * * * Intermidiate Functions * * * * * * * * * * * *//

void Foam::vofTwoPhaseCentralFoam::Initialize()
{
    const fvMesh & mesh = U_.mesh();
    const Foam::Time& runTime = U_.mesh().time();

    Info<< "\nReading g" << endl;
    const meshObjects::gravity& g = meshObjects::gravity::New(runTime);

    Info<< "\nReading hRef" << endl;
    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(dimLength, Zero)
    );

    Info<< "Calculating field g.h\n" << endl;
    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    gh_ = ((g & U_.mesh().C()) - ghRef);

    ghf_ = ((g & U_.mesh().Cf()) - ghRef);

    Compressibility();

    DensityThermo();

    MixtureProperties();

    updateLambda();

    speedOfSound();

    interpolateDensities();

    CharacteristicCourant();

    UpdateCentralWeights();
    Info << "UpdateCentralWeights();" << endl;
    UpdateCentralFields();
    Info << "UpdateCentralFields();" << endl;

    UEqn();
    Info << "UEqn();" << endl;

    // pressureGradient();

    dpdt_     = fvc::ddt(p_);
    TSource1_ = dpdt_;
    TSource2_ = dpdt_;

    p_rgh_.ref() = p_.internalField();
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryFieldRef(),
        phi_.boundaryField()
    );
    p_rghUpdatePatchFields();
    p_rgh_.correctBoundaryConditions();
    p_rgh_.write();
}

//* * * * * * * * * * * * * * * * * Flux` Functions * * * * * * * * * * * * *//

void Foam::vofTwoPhaseCentralFoam::updateKappa()
{
    bool kappaIsOne
    (
        pimple_.dict().getOrDefault("kappaIsOne", false)
    );

    bool kappaIsZero
    (
        pimple_.dict().getOrDefault("kappaIsZero", false)
    );

    if (kappaIsOne)
    {
        kappa_.primitiveFieldRef() = 1.0;
        kappa_.boundaryFieldRef() = 1.0;
    }

    if (kappaIsZero)
    {
        kappa_.primitiveFieldRef() = 0.0;
        kappa_.boundaryFieldRef() = 0.0;
    }

    if (!kappaIsOne && !kappaIsZero)
    {
        const fvMesh& mesh = U_.mesh();

        surfaceScalarField CfSf (max(CfSf_own_, CfSf_nei_));
        CfSf.setOriented(true);

        surfaceScalarField amaxSfbyDelta
        (
            mesh.surfaceInterpolation::deltaCoeffs()*amaxSf_
        );

        surfaceScalarField FaceAcCo
        (
            amaxSfbyDelta/mesh.magSf() * mesh.time().deltaT()
        );


        surfaceScalarField Maf (mag(phi_) / CfSf);

        kappa_ =
            min
            (
                Maf/FaceAcCo,
                scalar(1.0)
            );
    }

    onemkappa_ = 1.0 - kappa_;
    Info<< "max/min kappa: " << max(kappa_).value()
        << "/" << min(kappa_).value()
        << endl;

    //writeMaxMinKappa (kappa_);

    kappaBlend(kappa_, phi1_own_, phi1_nei_);
    kappaBlend(kappa_, phi2_own_, phi2_nei_);
    kappaBlend(kappa_, aphiv_own_, aphiv_nei_);
}

void Foam::vofTwoPhaseCentralFoam::writeMaxMinKappa
(
    const surfaceScalarField& kappa
)
{
    const fvMesh& mesh = kappa.mesh();
    if (kappa.mesh().time().outputTime())
    {
        volScalarField maxKappa
        (
            IOobject
            (
                "maxKappa",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("maxkappa", dimless,  0.0)
        );

        volScalarField minKappa
        (
            "minKappa",
            maxKappa*0.0
        );

        forAll(mesh.cells(), iCell)
        {
            scalar maxKappaCell =  -0.01;
            scalar minKappaCell =  1.01;

            const labelList& cellFaces = mesh.cells()[iCell];

            forAll(cellFaces, iFace)
            {
                if (mesh.isInternalFace(cellFaces[iFace]))
                {
                    if (kappa[cellFaces[iFace]] > maxKappaCell)
                    {
                        maxKappaCell = kappa[cellFaces[iFace]];
                    }
                    if (kappa[cellFaces[iFace]] < minKappaCell)
                    {
                        minKappaCell = kappa[cellFaces[iFace]];
                    }
                }
            }

            maxKappa[iCell] = maxKappaCell;
            minKappa[iCell] = minKappaCell;
        }

        maxKappa.write();
        minKappa.write();
    }
}

void Foam::vofTwoPhaseCentralFoam::combineMatrices
(
    const fvScalarMatrix& m1,
    const fvScalarMatrix& m2,
    const volScalarField& vf1,
    const volScalarField& vf2,
    fvScalarMatrix& m,
    bool removeConst
)
{
    autoPtr<lduMatrix> ldu1;
    autoPtr<lduMatrix> ldu2;

    if (removeConst)
    {
        ldu1.reset
        (
            new lduMatrix(const_cast<fvScalarMatrix&>(m1),true)
        );
        ldu2.reset
        (
            new lduMatrix(const_cast<fvScalarMatrix&>(m2),true)
        );
    }
    else
    {
        ldu1.reset(new lduMatrix(m1));
        ldu2.reset(new lduMatrix(m2));
    }

    ldu1().lduMatrix::operator*=
        (
            vf1.primitiveField()
        );
    ldu2().lduMatrix::operator*=
        (
            vf2.primitiveField()
        );
    m.lduMatrix::operator+=(ldu1());
    m.lduMatrix::operator+=(ldu2());

    m.source() +=
        vf1.primitiveField()*
        m1.source();
    m.source() +=
        vf2.primitiveField()*
        m2.source();

    forAll(m.boundaryCoeffs(), patchi)
    {
        scalarField pvf1
        (
             vf1.mesh().boundary()[patchi].patchInternalField(vf1.field())
        );
        scalarField pvf2
        (
             vf2.mesh().boundary()[patchi].patchInternalField(vf2.field())
        );

        m.internalCoeffs()[patchi] =
            pvf1*m1.internalCoeffs()[patchi] +
            pvf2*m2.internalCoeffs()[patchi];

        m.boundaryCoeffs()[patchi] =
            pvf1*m1.boundaryCoeffs()[patchi] +
            pvf2*m2.boundaryCoeffs()[patchi];
    }
}

void Foam::vofTwoPhaseCentralFoam::combineMatrices
(
    const fvVectorMatrix& m1,
    const fvVectorMatrix& m2,
    const volScalarField& vf1,
    const volScalarField& vf2,
    fvVectorMatrix& m,
    bool removeConst
)
{
    autoPtr<lduMatrix> ldu1;
    autoPtr<lduMatrix> ldu2;

    if (removeConst)
    {
        ldu1.reset
        (
            new lduMatrix(const_cast<fvVectorMatrix&>(m1),true)
        );
        ldu2.reset
        (
            new lduMatrix(const_cast<fvVectorMatrix&>(m2),true)
        );
    }
    else
    {
        ldu1.reset(new lduMatrix(m1));
        ldu2.reset(new lduMatrix(m2));
    }

    ldu1().lduMatrix::operator*=
        (
            vf1.primitiveField()
        );
    ldu2().lduMatrix::operator*=
        (
            vf2.primitiveField()
        );
    m.lduMatrix::operator+=(ldu1());
    m.lduMatrix::operator+=(ldu2());

    m.source() +=
        vf1.primitiveField()*
        m1.source();
    m.source() +=
        vf2.primitiveField()*
        m2.source();

    forAll(m.boundaryCoeffs(), patchi)
    {
        scalarField pvf1
        (
             vf1.mesh().boundary()[patchi].patchInternalField(vf1.field())
        );
        scalarField pvf2
        (
             vf2.mesh().boundary()[patchi].patchInternalField(vf2.field())
        );

        m.internalCoeffs()[patchi] =
            pvf1*m1.internalCoeffs()[patchi] +
            pvf2*m2.internalCoeffs()[patchi];

        m.boundaryCoeffs()[patchi] =
            pvf1*m1.boundaryCoeffs()[patchi] +
            pvf2*m2.boundaryCoeffs()[patchi];
    }
}

Foam::wordList Foam::vofTwoPhaseCentralFoam::p_rghPatchTypes()
{
    wordList patchTypes = p_.boundaryField().types();
    forAll(p_.boundaryField(), ipatch)
    {
        const fvPatchScalarField& pp =
            p_.boundaryField()[ipatch];
        if (pp.patchType() == Foam::fieldTypes::zeroGradientType) {
            patchTypes[ipatch] = fixedFluxPressureFvPatchScalarField::typeName;
        }
    }

    return patchTypes;
}

void Foam::vofTwoPhaseCentralFoam::p_rghUpdatePatchFields()
{
    forAll(p_.boundaryField(), ipatch)
    {
        const fvPatchScalarField& pp =
            p_.boundaryField()[ipatch];
        if (isA<totalPressureFvPatchScalarField>(pp))
        {
            totalPressureFvPatchScalarField& p_rghpf =
                dynamic_cast<totalPressureFvPatchScalarField&>
                (p_rgh_.boundaryFieldRef()[ipatch]);
            const totalPressureFvPatchScalarField& p_pf =
                dynamic_cast<const totalPressureFvPatchScalarField&>
                (p_.boundaryField()[ipatch]);
            p_rghpf.UName() = p_pf.UName();
            p_rghpf.phiName() = p_pf.phiName();
            p_rghpf.rhoName() = p_pf.rhoName();
            p_rghpf.psiName() = p_pf.psiName();
            p_rghpf.gamma() = p_pf.gamma();
            p_rghpf.p0() = p_pf.p0();
        }
    }
}

void Foam::vofTwoPhaseCentralFoam::DensityThermo()
{
    mixture_model_.density();
}

void Foam::vofTwoPhaseCentralFoam::MixtureProperties()
{
    mixture_model_.correct();
    turbulence_->correct();
}

void Foam::vofTwoPhaseCentralFoam::updateLambda()
{
    const auto& gam1 = mixture_model_.gamma1();
    const auto& gam2 = mixture_model_.gamma2();
    const auto& R1   = mixture_model_.R1();
    const auto& R2   = mixture_model_.R2();
    const auto& rho1 = mixture_model_.rho1();
    const auto& rho2 = mixture_model_.rho2();

    volScalarField C1sqr (gam1*R1*T_);
    volScalarField C2sqr (gam2*R2*T_);

    volScalarField Z1 (rho1*C1sqr);
    volScalarField Z2 (rho2*C2sqr);
    Lambda_ =
    (
        (volumeFraction1_*volumeFraction2_*(Z2 - Z1))
        /(Z1*volumeFraction2_ + Z2*volumeFraction1_)
    );
}

void Foam::vofTwoPhaseCentralFoam::Compressibility()
{
    mixture_model_.compressibility();
}

void Foam::vofTwoPhaseCentralFoam::speedOfSound()
{
    const auto& psi1 = mixture_model_.psi1();
    const auto& psi2 = mixture_model_.psi2();
    const auto& rho  = mixture_model_.rho();
    const auto& rho1 = mixture_model_.rho1();
    const auto& rho2 = mixture_model_.rho2();
    const auto& gam1 = mixture_model_.gamma1();
    const auto& gam2 = mixture_model_.gamma2();
    const auto& Cp1  = mixture_model_.Cp1();
    const auto& Cp2  = mixture_model_.Cp2();

    volScalarField psiM (volumeFraction1_*psi1 + volumeFraction2_*psi2);
    volScalarField y1 (volumeFraction1_*(rho1/rho));
    volScalarField y2 (volumeFraction2_*(rho2/rho));
    volScalarField CpM (y1*Cp1 + y2*Cp2);
    volScalarField CvM (y1*(Cp1/gam1) + y2*(Cp2/gam2));
    volScalarField gammaM (CpM/CvM);
    C_ = sqrt(gammaM/psiM);

    Cf_own_ = fvc::interpolate(C_, own_, "reconstruct(psi)");
    Cf_nei_ = fvc::interpolate(C_, nei_, "reconstruct(psi)");
}

void Foam::vofTwoPhaseCentralFoam::interpolateDensities()
{
    const auto& rho1 = mixture_model_.rho1();
    const auto& rho2 = mixture_model_.rho2();
    rho1_own_ =
        fvc::interpolate(rho1, own_, "reconstruct(" + rho1.name() + ")");
    rho1_nei_ =
        fvc::interpolate(rho1, nei_, "reconstruct(" + rho1.name() + ")");

    rho2_own_ =
        fvc::interpolate(rho2, own_, "reconstruct(" + rho2.name() + ")");
    rho2_nei_ =
        fvc::interpolate(rho2, nei_, "reconstruct(" + rho2.name() + ")");
}

void Foam::vofTwoPhaseCentralFoam::kappaBlend
(
    const surfaceScalarField& kappa,
    surfaceScalarField& flux_own,
    surfaceScalarField& flux_nei
)
{
    flux_own += (1.0 - kappa) * flux_nei;
    flux_nei = kappa * flux_nei;
}

