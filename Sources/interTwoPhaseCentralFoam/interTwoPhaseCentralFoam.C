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

#include "interTwoPhaseCentralFoam.H"

namespace Foam
{
    defineTypeNameAndDebug(interTwoPhaseCentralFoam, 0);
}

Foam::interTwoPhaseCentralFoam::interTwoPhaseCentralFoam(const fvMesh& mesh, pimpleControl& ctrl)
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

    R_
    (
        Foam::constant::physicoChemical::R
    ),

    R1_
    (
        dimensioned< scalar >("R1", *this)
    ),

    R2_
    (
        dimensioned< scalar >("R2", *this)
    ),

    molM1_
    (
        dimensioned< scalar >("molM1", *this)
    ),

    molM2_
    (
        dimensioned< scalar >("molM2", *this)
    ),

    psi1_
    (
        "psi1",
        1/(R1_*T_)
    ),

    psi2_
    (
        "psi2",
        1/(R2_*T_)
    ),

    rho1_
    (
        IOobject
        (
            "rho1",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
         ),
         psi1_*p_
    ),

    rho2_
    (
        IOobject
        (
            "rho2",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        psi2_*p_
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
        1 - volumeFraction1_
    ),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        volumeFraction1_*rho1_ + volumeFraction2_*rho2_
    ),

    rho01_
    (
        dimensioned< scalar >("rho01", *this)
    ),

    rho02_
    (
        dimensioned< scalar >("rho02", *this)
    ),

    Cp1_
    (
        dimensioned< scalar >("Cp1", *this)
    ),

    Cp2_
    (
        dimensioned< scalar >("Cp2", *this)
    ),

    mu1_
    (
        dimensioned< scalar >("mu1", *this)
    ),

    mu2_
    (
        dimensioned< scalar >("mu2", *this)
    ),

    Pr1_
    (
        dimensioned< scalar >("Pr1", *this)
    ),

    Pr2_
    (
        dimensioned< scalar >("Pr2", *this)
    ),

    alpha1_
    (
        "alpha1",
        mu1_/(Pr1_)
    ),

    alpha2_
    (
        "alpha2",
        mu2_/(Pr2_)
    ),

    gamma1_
    (
        dimensioned< scalar >("gamma1", Cp1_/(Cp1_ - (R_/molM1_)))
    ),

    gamma2_
    (
        dimensioned< scalar >("gamma2", Cp2_/(Cp2_ - (R_/molM2_)))
    ),

    C_
    (
        "C",
        sqrt(gamma1_*R1_*T_)
    ),

    K_
    (
        "K",
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

    rbyA_
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

/***********************Tadmor-Kurganov Scheme*******************************/

    v_zero
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

    alpha_own_
    (
        "alpha_own_ ",
        own_
    ),

    aSf_
    (
        "aSf_",
        (fvc::interpolate(U_, own_, "reconstruct(U)")) & mesh.Sf()
    ),

    alpha_nei_
    (
        "alpha_nei_",
        1.0 - alpha_own_
    ),

    alpha1_own_
    (
        "alpha1_own_ ",
        own_
    ),

    aSf1_
    (
        "aSf1_",
        (fvc::interpolate(U_, own_, "reconstruct(U)")) & mesh.Sf()
    ),

    alpha1_nei_
    (
        "alpha1_nei_",
        1.0 - alpha_own_
    ),

    alpha2_own_
    (
        "alpha2_own_ ",
        own_
    ),

    aSf2_
    (
        "aSf2_",
        (fvc::interpolate(U_, own_, "reconstruct(U)")) & mesh.Sf()
    ),

    alpha2_nei_
    (
        "alpha2_nei_",
        1.0 - alpha_own_
    ),

    phi1_own_
    (
        "phi1_own_",
        phi_*fvc::interpolate(rho1_, own_, "reconstruct(rho1)")*0.0
    ),

    phi1_nei_
    (
        "phi1_nei_",
        phi_*fvc::interpolate(rho1_, own_, "reconstruct(rho1)")*0.0
    ),

/*******************************Region Two**********************************/

    phi2_own_
    (
        "phi2_own_",
        phi_*fvc::interpolate(rho1_, own_, "reconstruct(rho1)")*0.0
    ),

    phi2_nei_
    (
        "phi2_nei_",
        phi_*fvc::interpolate(rho1_, own_, "reconstruct(rho1)")*0.0
    ),

/***********************Tadmor-Kurganov Scheme*******************************/

/*************************Pressure Equation*********************************/

    phi1d_own_
    (
        "phi1d_own_",
        phi_*fvc::interpolate(psi1_, own_, "reconstruct(psi)")
    ),

    phi1d_nei_
    (
        "phi1d_nei_",
        phi_*fvc::interpolate(psi1_, own_, "reconstruct(psi)")
    ),

    Dp1_own_
    (
        "Dp1_own_",
        alpha_own_ *fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)")
    ),

    Dp1_nei_
    (
        "Dp1_nei_",
        alpha_own_ *fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)")
    ),

    phi2d_own_
    (
        "phi2d_own_",
        phi_*fvc::interpolate(psi1_, own_, "reconstruct(psi)")
    ),

    phi2d_nei_
    (
        "phi2d_nei_",
        phi_*fvc::interpolate(psi1_, own_, "reconstruct(psi)")
    ),

    Dp2_own_
    (
        "Dp2_own_",
        alpha_own_ *fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)")
    ),

    Dp2_nei_
    (
        "Dp2_nei_",
        alpha_own_ *fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)")
    ),

    pEqn1_own_
    (
        fvm::div(phi1d_own_,p_) - fvm::laplacian(Dp1_own_, p_)
    ),

    pEqn1_nei_
    (
        fvm::div(phi1d_nei_,p_) - fvm::laplacian(Dp1_nei_, p_)
    ),

    pEqn2_own_
    (
        fvm::div(phi2d_own_,p_) - fvm::laplacian(Dp2_own_, p_)
    ),

    pEqn2_nei_
    (
        fvm::div(phi2d_nei_,p_) - fvm::laplacian(Dp2_nei_, p_)
    ),

    gradp_
    (
        "gradp_",
        fvc::grad(p_)
    ),

    divDevRhoReff1_
    (
        - fvm::laplacian(mu1_, U_)
        - fvc::div((mu1_)*dev2(Foam::T(fvc::grad(U_))))
    ),

    divDevRhoReff2_
    (
        - fvm::laplacian(mu2_, U_)
        - fvc::div((mu2_)*dev2(Foam::T(fvc::grad(U_))))
    ),

    Tviscosity1
    (
        - fvm::laplacian(alpha1_*Cp1_, T_)
    ),

    Tviscosity2
    (
        - fvm::laplacian(alpha2_*Cp2_, T_)
    ),

    devRhoReff1_
    (
        (-(alpha1_)*dev(twoSymm(fvc::grad(U_))))
    ),

    devRhoReff2_
    (
        (-(alpha2_)*dev(twoSymm(fvc::grad(U_))))
    ),

    TSource1_
    (
        fvc::ddt(p_)
    ),

    TSource2_
    (
        fvc::ddt(p_)
    ),

    TSource_
    (
        volumeFraction1_*1/Cp1_*TSource1_
    ),

    phi01d_own_
    (
        phi_*rho01_
    ),

    phi01d_nei_
    (
        phi_*rho01_
    ),

    phi02d_own_
    (
        phi_*rho02_
    ),

    phi02d_nei_
    (
        phi_*rho02_
    ),

    E1_
    (
        fvc::ddt(rho1_) + fvc::div(phi1_own_ + phi1_nei_)
    ),

    E2_
    (
        fvc::ddt(rho2_) + fvc::div(phi2_own_+phi2_nei_)
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

    rho1Min
    (
        dimensioned< scalar >("rho1Min", *this)
    ),

    rho2Min
    (
        dimensioned< scalar >("rho2Min", *this)
    ),

    Q_
    (
        0.5*magSqr(U_)
    )

{
    Info<< "\nConstructor is working\n" << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interTwoPhaseCentralFoam::~interTwoPhaseCentralFoam()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interTwoPhaseCentralFoam::saveOld()
{
    volumeFraction1_.oldTime();
    rho1_.oldTime();
    rho2_.oldTime();
    U_.oldTime();
    T_.oldTime();
    p_.oldTime();
    psi1_.oldTime();
    psi2_.oldTime();
    Q_.oldTime();
}


// * * * * * * * * * * * * * * * * Main Functions  * * * * * * * * * * * * * //


void Foam::interTwoPhaseCentralFoam::alpha1Eqnsolve()
{
    vF1face_ = fvc::interpolate(volumeFraction1_,"reconstruct(volumeFraction1)");
    vF2face_ = 1.0 - vF1face_;

    fvScalarMatrix alpha1Eqn
    (

        fvm::ddt(volumeFraction1_)
        +
        fvc::div(phi_, volumeFraction1_)
        ==
        (1 + K_)*fvc::div(phi_)*volumeFraction1_

    );

    alpha1Eqn.solve();

    volumeFraction2_ = 1 - volumeFraction1_;

    Info<< "max: volumeFraction1 " << max(volumeFraction1_).value()
        << " min: " << min(volumeFraction1_).value()
        << nl << endl;
}


void Foam::interTwoPhaseCentralFoam::UEqn()
{
    pressureGradient();
    Density();

    surfaceScalarField phiU_own = vF1face_*phi1_own_ + vF2face_*phi2_own_;
    surfaceScalarField phiU_nei = vF1face_*phi1_nei_ + vF2face_*phi2_nei_;

    E_ = fvc::ddt(rho_) + fvc::div(phiU_own) + fvc::div(phiU_nei);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_,U_) - fvm::Sp(E_,U_)
        +
        fvm::div(phiU_own,U_) + fvm::div(phiU_nei,U_)
    );

    rbyA_  = 1.0/UEqn.A();
    HbyA_ = UEqn.H()*rbyA_;
    HbyA_.boundaryFieldRef() == U_.boundaryField();
}

void Foam::interTwoPhaseCentralFoam::ReconstructVelocity()
{
    pressureGradient();
    U_ = HbyA_ - rbyA_*gradp_;
    U_.correctBoundaryConditions();
}


void Foam::interTwoPhaseCentralFoam::TEqnsolve()
{
    fvScalarMatrix TEqn
    (
        volumeFraction1_ *
        (
            fvm::ddt(rho1_,T_)
            + fvm::div(phi1_own_,T_) + fvm::div(phi1_nei_,T_)
//            - fvm::Sp(E1_,T_)
//          + Tviscosity1
            + 1/Cp1_*TSource1_
        )
        + volumeFraction2_*
        (
            fvm::ddt(rho2_,T_)
            + fvm::div(phi2_own_,T_) + fvm::div(phi2_nei_,T_)
//            - fvm::Sp(E2_,T_)
//          + Tviscosity2
            + 1/Cp2_*TSource2_
        )
    );

    TEqn.solve();
}


void Foam::interTwoPhaseCentralFoam::TEqnV2solve()
{
    surfaceScalarField phiU_own = vF1face_*phi1_own_ + vF2face_*phi2_own_;
    surfaceScalarField phiU_nei = vF1face_*phi1_nei_ + vF2face_*phi2_nei_;

    E_ = fvc::ddt(rho_) + fvc::div(phiU_own) + fvc::div(phiU_nei);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho_,T_)
        + fvm::div(phiU_own,T_) + fvm::div(phiU_nei,T_)
        - fvm::Sp(E_,T_)
        + TSource_
        + volumeFraction1_*1/Cp1_*TSource1_
        + volumeFraction2_*1/Cp2_*TSource2_
    );

    TEqn.solve();
}


void Foam::interTwoPhaseCentralFoam::pEqnsolve()
{
    pEqn1_own_ =
    (
        fvc::div(phi01d_own_)
        + fvm::div(phi1d_own_,p_)
        - fvm::laplacian(Dp1_own_, p_)
    );

    pEqn1_nei_ =
    (
        fvc::div(phi01d_nei_)
        + fvm::div(phi1d_nei_,p_)
        - fvm::laplacian(Dp1_nei_, p_)
    );

    pEqn2_own_ =
    (
        fvc::div(phi02d_own_)
        + fvm::div(phi2d_own_,p_)
        - fvm::laplacian(Dp2_own_, p_)
    );

    pEqn2_nei_ =
    (
        fvc::div(phi02d_nei_)
        + fvm::div(phi2d_nei_,p_)
        - fvm::laplacian(Dp2_nei_, p_)
    );

    fvScalarMatrix pEqn
    (
        volumeFraction1_*
        (
            fvm::ddt(psi1_,p_) + pEqn1_own_ + pEqn1_nei_
        )
        + volumeFraction2_*
        (
            fvm::ddt(psi2_,p_) + pEqn2_own_ + pEqn2_nei_
        )
    );

    pEqn.solve();
}


void Foam::interTwoPhaseCentralFoam::massError1()
{
    E1_ =
    (
        fvc::ddt(rho1_) + fvc::div(phi1_own_ + phi1_nei_)
    );

    Info<< " max: rho1 " << max(rho1_).value()
        << " min: " << min(rho1_).value()
        << nl << endl;
}


void Foam::interTwoPhaseCentralFoam::massError2()
{
    E2_ =
    (
        fvc::ddt(rho2_) + fvc::div(phi2_own_ + phi2_nei_)
    );

    Info<< " max: rho2 " << max(rho2_).value()
        << " min: " << min(rho2_).value()
        << nl << endl;
}


void Foam::interTwoPhaseCentralFoam::TSource()
{
    Q_ = 0.5*magSqr(U_);

    TSource1_ =
    (
      fvc::ddt(rho1_,Q_)
      + fvc::div(phi1_own_,Q_) + fvc::div(phi1_nei_,Q_)
      - fvc::ddt(p_)
//      - fvc::Sp(E1_,Q_)
    );

    TSource2_ =
    (
      fvc::ddt(rho2_,Q_)
      + fvc::div(phi2_own_,Q_) + fvc::div(phi2_nei_,Q_)
      - fvc::ddt(p_)
//      - fvc::Sp(E2_,Q_)
    );
}


void Foam::interTwoPhaseCentralFoam::TSourceV2()
{
    surfaceScalarField phiUCp_own = 1/Cp1_*vF1face_*phi1_own_ + 1/Cp2_*vF2face_*phi2_own_;
    surfaceScalarField phiUCp_nei = 1/Cp1_*vF1face_*phi1_nei_ + 1/Cp2_*vF2face_*phi2_nei_;

    Q_ = 0.5*magSqr(U_);

    TSource_ = fvc::div(phiUCp_own,Q_) + fvc::div(phiUCp_nei,Q_);

    TSource1_ =
    (
      fvc::ddt(rho1_,Q_)
      - fvc::ddt(p_)
      - fvc::Sp(E1_,Q_)
    );

    TSource2_ =
    (
      fvc::ddt(rho2_,Q_)
      - fvc::ddt(p_)
      - fvc::Sp(E2_,Q_)
    );
}

//* * * * * * * * * * * * * * * Intermidiate Functions * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::Initialize()
{
    R_ = Foam::constant::physicoChemical::R;

    Compressibility();

    DensityThermo();

    Density();

    updateK_();

    // Thermal conductivity
    alpha1_ = mu1_/Pr1_;
    alpha2_ = mu2_/Pr2_;

    UpdateCentralWeights();

//      divDevRhoReff();

    UEqn();

    UpdateCentralFields();

    pressureGradient();

//      devRhoReff();

//      Tviscosity1 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_));
//      Tviscosity2 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_));

    TSource1_ = 0*fvc::ddt(p_);
    TSource2_ = 0*fvc::ddt(p_);
    TSource_ = volumeFraction1_*1/Cp1_*TSource1_;
}

void Foam::interTwoPhaseCentralFoam::pressureGradient()
{
    surfaceScalarField p_own = fvc::interpolate(p_, own_, "reconstruct(p)");
    surfaceScalarField p_nei = fvc::interpolate(p_, nei_, "reconstruct(p)");

//    gradp_ = fvc::div((alpha_own_ *p_own + alpha_nei_*p_nei)*U_.mesh().Sf());
    gradp_ = fvc::div((alpha1_own_ *p_own + alpha1_nei_*p_nei)*U_.mesh().Sf());
//    gradp_ = fvc::grad(p_);
}

//* * * * * * * * * * * * * * * * * Flux` Functions * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::Flux()
{
    phi1_own_ = pEqn1_own_.flux() + phi01d_own_;
    phi1_nei_ = pEqn1_nei_.flux() + phi01d_nei_;

    phi2_own_ = pEqn2_own_.flux() + phi02d_own_;
    phi2_nei_ = pEqn2_nei_.flux() + phi02d_nei_;
}


void Foam::interTwoPhaseCentralFoam::volumeFlux()
{
  surfaceScalarField rho1_own = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
  surfaceScalarField rho1_nei = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

    phi_ =
    (
        (phi1_own_ + phi1_nei_)
        /(alpha_own_ *rho1_own + alpha_nei_*rho1_nei)
    );
}

//* * * * * * * * * * * * * * * Update Dencities * * * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::DensityThermo()
{
    rho1_ = psi1_*p_ + rho01_;
    rho2_ = psi2_*p_ + rho02_;
    rho1_ = max(rho1_,rho1Min);
    rho2_ = max(rho2_,rho2Min);
    rho1_.correctBoundaryConditions();
    rho2_.correctBoundaryConditions();
}


void Foam::interTwoPhaseCentralFoam::Density()
{
    rho_ = volumeFraction1_*rho1_ + volumeFraction2_*rho2_;
}

//* * * * * * * * * * * * * Update Dependable Variables * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::updateK_()
{
    volScalarField C1 = gamma1_*R1_*T_;
    volScalarField C2 = gamma2_*R2_*T_;

    volScalarField Z1 = rho1_*C1;
    volScalarField Z2 = rho2_*C2;

    K_ =
    (
        (volumeFraction1_*volumeFraction2_*(Z2 - Z1))
        /(Z1*volumeFraction2_ + Z2*volumeFraction1_)
    );
}


void Foam::interTwoPhaseCentralFoam::Compressibility()
{
    psi1_ = 1/(R1_*T_);
    psi2_ = 1/(R2_*T_);
}


void Foam::interTwoPhaseCentralFoam::speedOfSound()
{
    volScalarField rbypsiM =
        1/(volumeFraction1_*psi1_ + volumeFraction2_*psi2_);

    volScalarField psiM = 1/rbypsiM;

    volScalarField y1 = volumeFraction1_*(rho1_/rho_);

    volScalarField y2 = volumeFraction2_*(rho2_/rho_);

    volScalarField CpM = y1*Cp1_ + y2*Cp2_;

    volScalarField CvM = y1*(Cp1_/gamma1_) + y2*(Cp2_/gamma2_);

    volScalarField gammaM = CpM/CvM;

    C_ = sqrt(gammaM/psiM);
}

//* * * * * * * * * * * * Kurganov's coefficients Mixture * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::UpdateCentralWeights()
{

    surfaceScalarField rho1_own = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
    surfaceScalarField rho1_nei = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

    surfaceScalarField rho2_own = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
    surfaceScalarField rho2_nei = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

    surfaceScalarField phi1v_own = (fvc::interpolate(U_, own_, "reconstruct(U)")) & U_.mesh().Sf();
    surfaceScalarField phi1v_nei = (fvc::interpolate(U_, nei_, "reconstruct(U)")) & U_.mesh().Sf();

    speedOfSound();

    surfaceScalarField C_own = fvc::interpolate(C_, own_, "reconstruct(psi)");
    surfaceScalarField C_nei = fvc::interpolate(C_, nei_, "reconstruct(psi)");

    surfaceScalarField CSf_own = C_own*U_.mesh().magSf();
    CSf_own.setOriented(true);
    surfaceScalarField CSf_nei = C_nei*U_.mesh().magSf();
    CSf_nei.setOriented(true);

    surfaceScalarField ap = max(max(phi1v_own + CSf_own, phi1v_nei + CSf_nei), v_zero);
    surfaceScalarField am = min(min(phi1v_own - CSf_own, phi1v_nei - CSf_nei), v_zero);

    alpha_own_ = ap/(ap - am);
    aSf_ = am*alpha_own_ ;
    alpha_nei_ = 1.0 - alpha_own_ ;
    /*
    alpha_own_ = 0.5+0*ap/(ap - am);
    aSf_ = 0*am*alpha_own_;
    alpha_nei_ = 1.0 - alpha_own_ ;
    */
    }


void Foam::interTwoPhaseCentralFoam::UpdateCentralFields()
{
    surfaceScalarField rho1_own = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
    surfaceScalarField rho1_nei = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

    surfaceScalarField rho2_own = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
    surfaceScalarField rho2_nei = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

    surfaceScalarField psi1_own = fvc::interpolate(psi1_, own_, "reconstruct(psi1)");
    surfaceScalarField psi1_nei = fvc::interpolate(psi1_, nei_, "reconstruct(psi1)");

    surfaceScalarField psi2_own = fvc::interpolate(psi2_, own_, "reconstruct(psi2)");
    surfaceScalarField psi2_nei = fvc::interpolate(psi2_, nei_, "reconstruct(psi2)");

    surfaceVectorField rhoU1_own = fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)");
    surfaceVectorField rhoU1_nei = fvc::interpolate(rho1_*HbyA_, nei_, "reconstruct(U)");

    surfaceVectorField rhoU2_own = fvc::interpolate(rho2_*HbyA_, own_, "reconstruct(U)");
    surfaceVectorField rhoU2_nei = fvc::interpolate(rho2_*HbyA_, nei_, "reconstruct(U)");

    surfaceScalarField phi1v_own = (rhoU1_own/rho1_own) & U_.mesh().Sf();
    surfaceScalarField phi1v_nei = (rhoU1_nei/rho1_nei) & U_.mesh().Sf();

    surfaceScalarField phi2v_own = (rhoU2_own/rho2_own) & U_.mesh().Sf();
    surfaceScalarField phi2v_nei = (rhoU2_nei/rho2_nei) & U_.mesh().Sf();

    surfaceScalarField aphi1v_own = alpha_own_ *phi1v_own - aSf_;
    surfaceScalarField aphi1v_nei = alpha_nei_*phi1v_nei + aSf_;

    surfaceScalarField aphi2v_own = alpha_own_ *phi2v_own - aSf_;
    surfaceScalarField aphi2v_nei = alpha_nei_*phi2v_nei + aSf_;

    phi1d_own_ = aphi1v_own*psi1_own;
    phi1d_nei_ = aphi1v_nei*psi1_nei;

    phi2d_own_ = aphi2v_own*psi2_own;
    phi2d_nei_ = aphi2v_nei*psi2_nei;

    Dp1_own_ = alpha_own_ *fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)");
    Dp1_nei_ = alpha_nei_*fvc::interpolate(rho1_*rbyA_, nei_, "reconstruct(Dp)");

    Dp2_own_ = alpha_own_ *fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)");
    Dp2_nei_ = alpha_nei_*fvc::interpolate(rho2_*rbyA_, nei_, "reconstruct(Dp)");

    phi01d_own_ = aphi1v_own*rho01_;
    phi01d_nei_ = aphi1v_nei*rho01_;

    phi02d_own_ = aphi2v_own*rho02_;
    phi02d_nei_ = aphi2v_nei*rho02_;
}

//* * * * * * * * * Kurganov's coefficients Individual Phases * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::UpdateCentralWeightsIndividual()
{
    surfaceScalarField phi1v_own = (fvc::interpolate(U_, own_, "reconstruct(U)")) & U_.mesh().Sf();
    surfaceScalarField phi1v_nei = (fvc::interpolate(U_, nei_, "reconstruct(U)")) & U_.mesh().Sf();

//************************************ Phase One *****************************//

    surfaceScalarField rho1_own = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
    surfaceScalarField rho1_nei = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

    volScalarField C1 = sqrt(gamma1_*R1_*T_);

    surfaceScalarField C1_own = fvc::interpolate(C1, own_, "reconstruct(psi)");
    surfaceScalarField C1_nei = fvc::interpolate(C1, nei_, "reconstruct(psi)");

    surfaceScalarField C1Sf_own = C1_own*U_.mesh().magSf();
    C1Sf_own.setOriented(true);
    surfaceScalarField C1Sf_nei = C1_nei*U_.mesh().magSf();
    C1Sf_nei.setOriented(true);

    surfaceScalarField ap1 = max(max(phi1v_own + C1Sf_own, phi1v_nei + C1Sf_nei), v_zero);
    surfaceScalarField am1 = min(min(phi1v_own - C1Sf_own, phi1v_nei - C1Sf_nei), v_zero);

    alpha1_own_ = ap1/(ap1 - am1);
    aSf1_ = am1*alpha1_own_ ;
    alpha1_nei_ = 1.0 - alpha1_own_ ;

    /*
    alpha1_own_ = 0.5+0*ap1/(ap1 - am1);
    aSf1_ = 0*am1*alpha1_own_ ;
    alpha1_nei_ = 1.0 - alpha1_own_ ;
    */
//************************************ Phase Two *****************************//

    surfaceScalarField rho2_own = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
    surfaceScalarField rho2_nei = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

    volScalarField C2 = sqrt(gamma2_*R2_*T_);

    surfaceScalarField C2_own = fvc::interpolate(C2, own_, "reconstruct(psi)");
    surfaceScalarField C2_nei = fvc::interpolate(C2, nei_, "reconstruct(psi)");

    surfaceScalarField C2Sf_own = C2_own*U_.mesh().magSf();
    C2Sf_own.setOriented(true);
    surfaceScalarField C2Sf_nei = C2_nei*U_.mesh().magSf();
    C2Sf_nei.setOriented(true);

    surfaceScalarField ap2 = max(max(phi1v_own + C2Sf_own, phi1v_nei + C2Sf_nei), v_zero);
    surfaceScalarField am2 = min(min(phi1v_own - C2Sf_own, phi1v_nei - C2Sf_nei), v_zero);

    alpha2_own_ = ap2/(ap2 - am2);
    aSf2_ = am2*alpha2_own_ ;
    alpha2_nei_ = 1.0 - alpha2_own_ ;

    Info<< " C1 " << max(C1).value()
        << " C2 " << max(C2).value()
        << nl << endl;
    }


void Foam::interTwoPhaseCentralFoam::UpdateCentralFieldsIndividual()
{
    surfaceScalarField rho1_own = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
    surfaceScalarField rho1_nei = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

    surfaceScalarField rho2_own = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
    surfaceScalarField rho2_nei = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

    surfaceScalarField psi1_own = fvc::interpolate(psi1_, own_, "reconstruct(psi1)");
    surfaceScalarField psi1_nei = fvc::interpolate(psi1_, nei_, "reconstruct(psi1)");

    surfaceScalarField psi2_own = fvc::interpolate(psi2_, own_, "reconstruct(psi2)");
    surfaceScalarField psi2_nei = fvc::interpolate(psi2_, nei_, "reconstruct(psi2)");

    surfaceVectorField rhoU1_own = fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)");
    surfaceVectorField rhoU1_nei = fvc::interpolate(rho1_*HbyA_, nei_, "reconstruct(U)");

    surfaceVectorField rhoU2_own = fvc::interpolate(rho2_*HbyA_, own_, "reconstruct(U)");
    surfaceVectorField rhoU2_nei = fvc::interpolate(rho2_*HbyA_, nei_, "reconstruct(U)");

    surfaceScalarField phi1v_own = (rhoU1_own/rho1_own) & U_.mesh().Sf();
    surfaceScalarField phi1v_nei = (rhoU1_nei/rho1_nei) & U_.mesh().Sf();

    surfaceScalarField phi2v_own = (rhoU2_own/rho2_own) & U_.mesh().Sf();
    surfaceScalarField phi2v_nei = (rhoU2_nei/rho2_nei) & U_.mesh().Sf();

    surfaceScalarField aphi1v_own = alpha1_own_ *phi1v_own - aSf1_;
    surfaceScalarField aphi1v_nei = alpha1_nei_*phi1v_nei + aSf1_;

    surfaceScalarField aphi2v_own = alpha2_own_ *phi2v_own - aSf2_;
    surfaceScalarField aphi2v_nei = alpha2_nei_*phi2v_nei + aSf2_;

    phi1d_own_ = aphi1v_own*psi1_own;
    phi1d_nei_ = aphi1v_nei*psi1_nei;

    phi2d_own_ = aphi2v_own*psi2_own;
    phi2d_nei_ = aphi2v_nei*psi2_nei;

    Dp1_own_ = alpha1_own_ *fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)");
    Dp1_nei_ = alpha1_nei_*fvc::interpolate(rho1_*rbyA_, nei_, "reconstruct(Dp)");

    Dp2_own_ = alpha2_own_ *fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)");
    Dp2_nei_ = alpha2_nei_*fvc::interpolate(rho2_*rbyA_, nei_, "reconstruct(Dp)");

    phi01d_own_ = aphi1v_own*rho01_;
    phi01d_nei_ = aphi1v_nei*rho01_;

    phi02d_own_ = aphi2v_own*rho02_;
    phi02d_nei_ = aphi2v_nei*rho02_;
}

//* * * * * * * * * * * * * * * * * Viscosity * * * * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::divDevRhoReff()
{
    divDevRhoReff1_ =
    (
        - fvm::laplacian(mu1_, U_)
        - fvc::div((mu1_)*dev2(Foam::T(fvc::grad(U_))))
    );

    divDevRhoReff2_ =
    (
        - fvm::laplacian(mu2_, U_)
        - fvc::div((mu2_)*dev2(Foam::T(fvc::grad(U_))))
    );
}


void Foam::interTwoPhaseCentralFoam::viscosityTEqn()
{
    Tviscosity1 = - fvm::laplacian(alpha1_*Cp1_, T_);

    Tviscosity2 =- fvm::laplacian(alpha2_*Cp2_, T_);
}


void Foam::interTwoPhaseCentralFoam::devRhoReff()
{
    devRhoReff1_ = (-(alpha1_)*dev(twoSymm(fvc::grad(U_))));
    devRhoReff2_ = (-(alpha2_)*dev(twoSymm(fvc::grad(U_))));
    // Rename TSourse
    TSource1_ =
    fvc::div((linearInterpolate((-devRhoReff1_) & U_) & U_.mesh().Sf())());

    TSource2_ =
    fvc::div((linearInterpolate((-devRhoReff2_) & U_) & U_.mesh().Sf())());
}

//* * * * * * * * * * * * * * * * * * Others * * * * * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::divU()
{

    surfaceScalarField rbyAf = fvc::interpolate(rbyA_);
/*
    Foam::CorrectPhi
    (
        U_,
        phi_,
        p_,
        rbyAf,
        divU_,
        pimple_
    );
*/
}
