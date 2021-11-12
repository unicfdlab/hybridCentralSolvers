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
      dimensionedScalar
      (
        "R",
        dimensionSet(1, 2, -2, -1, -1, 0, 0),
        8.31446261815324
      )
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

    rho01Sc_
    (
      dimensioned< scalar >("rho01", *this)
    ),

    rho02Sc_
    (
      dimensioned< scalar >("rho02", *this)
    ),

    rho01_
    (
      IOobject
      (
          "rho2",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      (rho01Sc_*rho1_)/rho1_
    ),

    rho02_
    (
      IOobject
      (
          "rho2",
          mesh.time().timeName(),
          mesh,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      (rho02Sc_*rho2_)/rho2_
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
      dimensioned< scalar >("gamma1", Cp1_/(Cp1_-(R_/molM1_)))
    ),

    gamma2_
    (
      dimensioned< scalar >("gamma2", Cp2_/(Cp2_-(R_/molM2_)))
    ),

    C1_
    (
      "C1",
      gamma1_*R1_*T_
    ),

    C2_
    (
      "C2",
      gamma2_*R2_*T_
    ),

    C_
    (
      "C",
      volumeFraction1_*sqrt(C1_) + volumeFraction2_*sqrt(C2_)
    ),

    Z1_
    (
      "Z1",
      rho1_*C1_
    ),

    Z2_
    (
      "Z2",
      rho2_*C2_
    ),

    K_
    (
      "K",
      (volumeFraction1_*volumeFraction2_*(Z1_ - Z2_))
      /
      (Z1_*volumeFraction2_ + Z2_*volumeFraction1_)
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

    rho1AmpU_
    (
      "rho1AmplitudeU",
      magSqr(U_)*rho1_
    ),

    rho2AmpU_
    (
      "rho2AmplitudeU",
      magSqr(U_)*rho2_
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

    rho1_own
    (
        "rho1own",
        fvc::interpolate(rho1_, own_, "reconstruct(rho1)")
    ),

    rho1_nei
    (
        "rho1nei",
        fvc::interpolate(rho1_, nei_, "reconstruct(rho1)")
    ),

    phi1v_own
    (
        "phi1v_own",
        (fvc::interpolate(rho1_*U_, own_, "reconstruct(U)")/rho1_own) & mesh.Sf()
    ),

    phi1v_nei
    (
        "phi1v_nei",
        (fvc::interpolate(rho1_*U_, nei_, "reconstruct(U)")/rho1_nei) & mesh.Sf()
    ),

     C_own
    (
        "C_own",
        fvc::interpolate(sqrt(C1_), own_, "reconstruct(phi)")
    ),

     C_nei
    (
        "C_nei",
        fvc::interpolate(sqrt(C1_), nei_, "reconstruct(phi)")
    ),

     CSf_own
    (
        "CSf_own",
        C_own * mesh.magSf()
    ),

     CSf_nei
    (
        "CSf_nei",
        C_nei * mesh.magSf()
    ),

     ap
    (
      IOobject
      (
          "ap",
          mesh.time().timeName(),
          mesh
      ),
      mesh,
      dimensionedScalar("ap", dimensionSet(0, 3, -1, 0, 0, 0, 0), 1.0)
    ),

     am
    (
      IOobject
      (
          "am",
          mesh.time().timeName(),
          mesh
      ),
      mesh,
      dimensionedScalar("ap", dimensionSet(0, 3, -1, 0, 0, 0, 0), -1.0)
    ),

     alpha_own
    (
        "alpha_own ",
        ap/(ap - am)
    ),

     aSf
    (
        "aSf",
        am*alpha_own
    ),

     alpha_nei
    (
        "alpha_nei",
        1.0 - alpha_own
    ),

     phi1_own
    (
        "phi1_own",
        phi_*rho1_own * 0.0
    ),

     phi1_nei
    (
        "phi1_nei",
        phi_*rho1_own * 0.0
    ),

/*******************************Region Two**********************************/

    rho2_own
    (
        "rho2own",
        fvc::interpolate(rho2_, own_, "reconstruct(rho2)")
    ),

    rho2_nei
    (
        "rho2nei",
        fvc::interpolate(rho2_, nei_, "reconstruct(rho2)")
    ),

    phi2v_own
    (
        "phi2v_own",
        (fvc::interpolate(rho1_*U_, own_, "reconstruct(U)")/rho1_own) & mesh.Sf()
    ),

    phi2v_nei
    (
        "phi2v_nei",
        (fvc::interpolate(rho1_*U_, nei_, "reconstruct(U)")/rho1_nei) & mesh.Sf()
    ),

     phi2_own
    (
        "phi2_own",
        phi_*rho2_own * 0.0
    ),

     phi2_nei
    (
        "phi2_nei",
        phi_*rho2_own * 0.0
    ),

/***********************Tadmor-Kurganov Scheme*******************************/

/*************************Pressure Equation*********************************/

    psi1_own
    (
        "psi1_own",
        fvc::interpolate(psi1_, own_, "reconstruct(psi)")
    ),

    psi1_nei
    (
        "psi1_nei",
        fvc::interpolate(psi1_, nei_, "reconstruct(psi)")
    ),

    psiU1_own
    (
        "psiU1_own",
         fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)")

    ),

    psiU1_nei
    (
        "psiU1_nei",
        fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)")
    ),

    aphi1v_own
    (
      "aphi1v_own",
      phi1v_own
    ),

    aphi1v_nei
    (
      "aphi1v_nei",
      phi1v_nei
    ),

    phi1d_own
    (
        "phi1d_own",
        aphi1v_own * psi1_own
    ),

    phi1d_nei
    (
        "phi1d_nei",
        aphi1v_nei * psi1_nei
    ),

    Dp1_own
    (
        "Dp1_own",
        alpha_own *fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)")
    ),

    Dp1_nei
    (
        "Dp1_nei",
        alpha_own *fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)")
    ),

    psi2_own
    (
        "psi2_own",
        fvc::interpolate(psi2_, own_, "reconstruct(psi)")
    ),

    psi2_nei
    (
        "psi2_nei",
        fvc::interpolate(psi2_, nei_, "reconstruct(psi)")
    ),

    psiU2_own
    (
        "psiU2_own",
        fvc::interpolate(rho2_*HbyA_, own_, "reconstruct(U)")
    ),

    psiU2_nei
    (
        "psiU2_nei",
         fvc::interpolate(rho2_*HbyA_, nei_, "reconstruct(U)")
    ),

    aphi2v_own
    (
      "aphi2v_own",
      phi1v_own
    ),

    aphi2v_nei
    (
      "aphi2v_nei",
      phi1v_nei
    ),

    phi2d_own
    (
        "phi2d_own",
        aphi2v_own * psi2_own
    ),

    phi2d_nei
    (
        "phi2d_nei",
        aphi2v_nei * psi2_nei
    ),

    Dp2_own
    (
        "Dp2_own",
        alpha_own *fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)")
    ),

    Dp2_nei
    (
        "Dp2_nei",
        alpha_own *fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)")
    ),

    pEqn1_own
    (
        fvm::div(phi1d_own,p_) - fvm::laplacian(Dp1_own, p_)
    ),

    pEqn1_nei
    (
        fvm::div(phi1d_nei,p_) - fvm::laplacian(Dp1_nei, p_)
    ),

    pEqn2_own
    (
        fvm::div(phi2d_own,p_) - fvm::laplacian(Dp2_own, p_)
    ),

    pEqn2_nei
    (
        fvm::div(phi2d_nei,p_) - fvm::laplacian(Dp2_nei, p_)
    ),

    p_own
    (
        "p_own",
        fvc::interpolate(p_, own_, "reconstruct(p)")
    ),

    p_nei
    (
        "p_nei",
        fvc::interpolate(p_, own_, "reconstruct(p)")
    ),

    gradp_
    (
      "gradp_",
      fvc::div((alpha_own *p_own + alpha_nei*p_nei)*mesh.Sf())
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
       - 0.5*fvc::laplacian(alpha1_, magSqr(U_))
     ),

     Tviscosity2
     (
       - fvm::laplacian(alpha2_*Cp2_, T_)
       - 0.5*fvc::laplacian(alpha2_, magSqr(U_))
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
       fvc::ddt(p_) - 0.5*fvc::ddt(rho1AmpU_)
     ),

     TSource2_
     (
       fvc::ddt(p_) - 0.5*fvc::ddt(rho2AmpU_)
     ),

     psi1f
     (
       "psi1f",
       linearInterpolate(psi1_)
     ),

     psi2f
     (
       "psi2f",
       linearInterpolate(psi1_)
     ),

     phi01d_own
     (
       aphi1v_own*rho01Sc_
     ),

     phi01d_nei
     (
       aphi1v_nei*rho01Sc_
     ),

     phi02d_own
     (
       aphi2v_own*rho02Sc_
     ),

     phi02d_nei
     (
       aphi2v_nei*rho02Sc_
     ),

     E1_
     (
        fvc::ddt(rho1_) + fvc::div(phi1_own + phi1_nei)
     ),

     E2_
     (
         fvc::ddt(rho2_) + fvc::div(phi2_own+phi2_nei)
     ),

     E_
     (
        volumeFraction1_*E1_ + volumeFraction2_*E2_
     ),

     ddtvF1_
     (
       fvc::ddt(volumeFraction1_)
     ),

     ddtvF2_
     (
       fvc::ddt(volumeFraction2_)
     ),

     vF1face_
     (
       fvc::interpolate(volumeFraction1_, "reconstruct(volumeFraction1)")
     ),

     vF2face_
     (
       1.0 - vF1face_
     ),

     vFPhi1_
     (
       phi_*vF1face_
     ),

     vFPhi2_
     (
       phi_ - vFPhi1_
     ),

     phiU_own
     (
       "phiU_own",
        vF1face_*phi1_own + vF2face_*phi2_own
     ),

     phiU_nei
     (
       "phiU_nei",
        vF1face_*phi1_nei + vF2face_*phi2_nei
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
      rho01_.oldTime();
      rho02_.oldTime();
      rho1AmpU_.oldTime();
      rho2AmpU_.oldTime();
      TSource1_.oldTime();
      TSource2_.oldTime();
      ddtvF1_.oldTime();
      ddtvF2_.oldTime();

    }


// * * * * * * * * * * * * * * * * Main Functions  * * * * * * * * * * * * * //


    void Foam::interTwoPhaseCentralFoam::alpha1Eqnsolve()
    {

      vF1face_ = fvc::interpolate(volumeFraction1_,
          "reconstruct(volumeFraction1)");
      vF2face_ = 1.0 - vF1face_;
      vFPhi1_ = phi_*vF1face_;
      vFPhi2_ = phi_ - vFPhi1_;

      fvScalarMatrix alpha1Eqn
      (

        fvm::ddt(volumeFraction1_)
        +
        fvc::div(phi_, volumeFraction1_)
        ==
        (1+K_)*fvc::div(phi_)*volumeFraction1_

      );

      alpha1Eqn.solve();

      volumeFraction2_ = 1 - volumeFraction1_;


      ddtvF1_ = fvc::ddt(volumeFraction1_);
      ddtvF2_ = fvc::ddt(volumeFraction2_);

      Info<< "max: volumeFraction1 " << max(volumeFraction1_).value()
          << " min: " << min(volumeFraction1_).value()
          << nl << endl;
    }


    void Foam::interTwoPhaseCentralFoam::UEqn()
    {

      pressureGradient();
      Density();

      phiU_own = vF1face_*phi1_own+vF2face_*phi2_own;
      phiU_nei = vF1face_*phi1_nei+vF2face_*phi2_nei;

      E_ = fvc::ddt(rho_) + fvc::div(phiU_own) + fvc::div(phiU_nei);

      fvVectorMatrix UEqn
      (
        fvm::ddt(rho_,U_) - fvm::Sp(E_,U_)
        +
        fvm::div(phiU_own,U_) + fvm::div(phiU_nei,U_)
      );

      rbyA_  = 1.0 / UEqn.A();
      HbyA_ = UEqn.H() * rbyA_;
      HbyA_.boundaryFieldRef() == U_.boundaryField();
    }


    void Foam::interTwoPhaseCentralFoam::TEqnsolve()
    {

      fvScalarMatrix TEqn
      (
        volumeFraction1_ *
        (
          Cp1_*fvm::ddt(rho1_,T_)
          -
          Cp1_*fvm::Sp(E1_,T_)
          +
          Cp1_*fvm::div(phi1_own,T_) + Cp1_*fvm::div(phi1_nei,T_)
          +
          0.5*(fvc::div(phi1_own,magSqr(U_)) + fvc::div(phi1_nei,magSqr(U_)))
//          +
//          Tviscosity1
          -
          TSource1_
        )
        +
        volumeFraction2_ *
        (
          Cp2_*fvm::ddt(rho2_,T_)
          -
          Cp2_*fvm::Sp(E2_,T_)
          +
          Cp2_*fvm::div(phi2_own,T_) + Cp2_*fvm::div(phi2_nei,T_)
          +
          0.5*(fvc::div(phi2_own,magSqr(U_)) + fvc::div(phi2_nei,magSqr(U_)))
//          +
//          Tviscosity2
          -
          TSource2_
        )
      );

      TEqn.solve();

    }


    void Foam::interTwoPhaseCentralFoam::rho1Eqnsolve()
    {

      E1_ =
      (
        fvc::ddt(rho01_)+fvc::ddt(psi1_,p_)
        +
        fvc::div(phi1_own+phi1_nei)
      );

      Info<< " max: rho1 " << max(rho1_).value()
          << " min: " << min(rho1_).value()
          << nl << endl;
    }


    void Foam::interTwoPhaseCentralFoam::rho2Eqnsolve()
    {

      E2_ =
      (
        fvc::ddt(rho02_)+fvc::ddt(psi2_,p_)
        +
        fvc::div(phi2_own+phi2_nei)
      );

      Info<< " max: rho2 " << max(rho2_).value()
          << " min: " << min(rho2_).value()
          << nl << endl;
    }


    void Foam::interTwoPhaseCentralFoam::pEqnsolve()
    {

      pEqn1_own =
      (
          fvc::div(phi01d_own)
          + fvm::div(phi1d_own,p_)
          - fvm::laplacian(Dp1_own, p_)
      );

      pEqn1_nei =
      (
          fvc::div(phi01d_nei)
          + fvm::div(phi1d_nei,p_)
          - fvm::laplacian(Dp1_nei, p_)
      );

      pEqn2_own =
      (
          fvc::div(phi02d_own)
          + fvm::div(phi2d_own,p_)
          - fvm::laplacian(Dp2_own, p_)
      );

      pEqn2_nei =
      (
          fvc::div(phi02d_nei)
          + fvm::div(phi2d_nei,p_)
          - fvm::laplacian(Dp2_nei, p_)
      );

      fvScalarMatrix pEqn
      (
          volumeFraction1_*
          (
            fvm::ddt(psi1_,p_) + fvc::ddt(rho01_)
            + pEqn1_own + pEqn1_nei
          )
          +
          volumeFraction2_*
          (
            fvm::ddt(psi2_,p_) + fvc::ddt(rho02_)
            + pEqn2_own + pEqn2_nei
          )
      );

      pEqn.solve();

    }

//* * * * * * * * * * * * * * * Intermidiate Functions * * * * * * * * * * * *//

     void Foam::interTwoPhaseCentralFoam::Flux()
     {

        phi1_own = pEqn1_own.flux() + phi01d_own;
        phi1_nei = pEqn1_nei.flux() + phi01d_nei;

        phi2_own = pEqn2_own.flux() + phi02d_own;
        phi2_nei = pEqn2_nei.flux() + phi02d_nei;

      }



      void Foam::interTwoPhaseCentralFoam::CompressibilityCoefficient()
      {

        C1_ = gamma1_*R1_*T_;
        C2_ = gamma2_*R2_*T_;

        Z1_ = rho1_*C1_;
        Z2_ = rho2_*C2_;

        K_ = (volumeFraction1_*volumeFraction2_*(Z2_ - Z1_))
        /(Z1_*volumeFraction2_ + Z2_*volumeFraction1_);

      }


    void Foam::interTwoPhaseCentralFoam::volumeFlux()
    {

      phi_ =
      (phi1_own + phi1_nei)/(alpha_own *rho1_own + alpha_nei*rho1_nei);

      Info << "C1: " << max(sqrt(C1_)).value() << endl;
      Info << "C2: " << max(sqrt(C2_)).value() << endl;
      Info << "C: " << max(sqrt(C_)).value() << endl;

    }


    void Foam::interTwoPhaseCentralFoam::ReconstructVelocity()
    {

      pressureGradient();
      U_ = HbyA_ - rbyA_*gradp_;
      U_.correctBoundaryConditions();

    }

    void Foam::interTwoPhaseCentralFoam::DensityThermo()
    {

      dimensionedScalar rho2Min ("rho2Min", rho2_.dimensions(), 0.01 );
      rho1_ = psi1_*p_ + rho01Sc_;
      rho2_ = psi2_*p_ + rho02Sc_;
      rho2_ = max(rho2_,rho2Min);
      rho1_.correctBoundaryConditions();
      rho2_.correctBoundaryConditions();

    }


    void Foam::interTwoPhaseCentralFoam::Density()
    {

      rho_ = volumeFraction1_*rho1_ + volumeFraction2_*rho2_;

    }

    void Foam::interTwoPhaseCentralFoam::Compressibility()
    {

      psi1_ = 1/(R1_*T_);
      psi2_ = 1/(R2_*T_);

    }

    void Foam::interTwoPhaseCentralFoam::Initialize()
    {

      psi1_ = (1/(R1_*T_));
      psi2_ = (1/(R2_*T_));

      rho01_ = (rho01Sc_*rho1_)/rho1_;
      rho02_ = (rho02Sc_*rho2_)/rho2_;

      rho1_ = psi1_*p_ + rho01_;
      rho2_ = psi2_*p_ + rho02_;
      rho1_.correctBoundaryConditions();
      rho2_.correctBoundaryConditions();

      rho_ = volumeFraction1_*rho1_ + volumeFraction2_*rho2_;

      gamma1_ = Cp1_/(Cp1_ - (R_/molM1_));
      gamma2_ = Cp2_/(Cp2_ - (R_/molM2_));

      CompressibilityCoefficient();

      ThermalConductivity();
      UpdateCentralWeights();

//      divDevRhoReff();

      UEqn();

      UpdateCentralFields();

      pressureGradient();

//      devRhoReff();

      rho1AmpU_ = magSqr(U_)*rho1_;
      rho2AmpU_ = magSqr(U_)*rho2_;

//      Tviscosity1 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_));
//      Tviscosity2 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_));

      TSource1_ = 0*(fvc::ddt(p_) - 0.5*fvc::ddt(rho1AmpU_));

      TSource2_ = 0*(fvc::ddt(p_) - 0.5*fvc::ddt(rho2AmpU_));

    }


    void Foam::interTwoPhaseCentralFoam::pressureGradient()
    {

      p_own = fvc::interpolate(p_, own_, "reconstruct(p)");
      p_nei = fvc::interpolate(p_, nei_, "reconstruct(p)");

      gradp_ = fvc::div((alpha_own *p_own + alpha_nei*p_nei)*U_.mesh().Sf());

    }


    void Foam::interTwoPhaseCentralFoam::UpdateCentralWeights()
    {

      rho1_own = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
      rho1_nei = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

      rho2_own = fvc::interpolate(rho2_, own_, "reconstruct(rho2)");
      rho2_nei = fvc::interpolate(rho2_, nei_, "reconstruct(rho2)");

      phi1v_own = (fvc::interpolate(U_, own_, "reconstruct(U)")) & U_.mesh().Sf();
      phi1v_nei = (fvc::interpolate(U_, nei_, "reconstruct(U)")) & U_.mesh().Sf();

      speedOfSound();

      C_own      = fvc::interpolate((C_), own_, "reconstruct(psi1)");
      C_nei      = fvc::interpolate((C_), nei_, "reconstruct(psi1)");


      CSf_own = C_own*U_.mesh().magSf();
      CSf_own.setOriented(true);
      CSf_nei = C_nei*U_.mesh().magSf();
      CSf_nei.setOriented(true);


      ap = max(max(phi1v_own + CSf_own, phi1v_nei + CSf_nei), v_zero);
      am = min(min(phi1v_own - CSf_own, phi1v_nei - CSf_nei), v_zero);

      alpha_own    = ap/(ap - am);
      aSf         = am*alpha_own ;
      alpha_nei   = 1.0 - alpha_own ;

    }


    void Foam::interTwoPhaseCentralFoam::UpdateCentralFields()
    {

      rho1_own = fvc::interpolate(rho1_, own_, "reconstruct(rho1)");
      rho1_nei = fvc::interpolate(rho1_, nei_, "reconstruct(rho1)");

      rho2_own = fvc::interpolate(rho2_, own_, "reconstruct(rho1)");
      rho2_nei = fvc::interpolate(rho2_, nei_, "reconstruct(rho1)");

      psi1_own = fvc::interpolate(psi1_, own_, "reconstruct(psi)");
      psi1_nei = fvc::interpolate(psi1_, nei_, "reconstruct(psi)");

      psi2_own = fvc::interpolate(psi2_, own_, "reconstruct(psi)");
      psi2_nei = fvc::interpolate(psi2_, nei_, "reconstruct(psi)");

      psiU1_own = fvc::interpolate(rho1_*HbyA_, own_, "reconstruct(U)");
      psiU1_nei = fvc::interpolate(rho1_*HbyA_, nei_, "reconstruct(U)");

      psiU2_own = fvc::interpolate(rho2_*HbyA_, own_, "reconstruct(U)");
      psiU2_nei = fvc::interpolate(rho2_*HbyA_, nei_, "reconstruct(U)");

      phi1v_own = (psiU1_own/rho1_own) & U_.mesh().Sf();
      phi1v_nei = (psiU1_nei/rho1_nei) & U_.mesh().Sf();

      phi2v_own = (psiU2_own/rho2_own) & U_.mesh().Sf();
      phi2v_nei = (psiU2_nei/rho2_nei) & U_.mesh().Sf();

      aphi1v_own = alpha_own *phi1v_own - aSf;
      aphi1v_nei = alpha_nei*phi1v_nei + aSf;

      aphi2v_own = alpha_own *phi2v_own - aSf;
      aphi2v_nei = alpha_nei*phi2v_nei + aSf;

      phi1d_own = aphi1v_own*psi1_own;
      phi1d_nei = aphi1v_nei*psi1_nei;

      phi2d_own = aphi2v_own*psi2_own;
      phi2d_nei = aphi2v_nei*psi2_nei;

      Dp1_own = alpha_own *fvc::interpolate(rho1_*rbyA_, own_, "reconstruct(Dp)");
      Dp1_nei = alpha_nei*fvc::interpolate(rho1_*rbyA_, nei_, "reconstruct(Dp)");

      Dp2_own = alpha_own *fvc::interpolate(rho2_*rbyA_, own_, "reconstruct(Dp)");
      Dp2_nei = alpha_nei*fvc::interpolate(rho2_*rbyA_, nei_, "reconstruct(Dp)");

      phi01d_own = aphi1v_own*rho01Sc_;
      phi01d_nei = aphi1v_nei*rho01Sc_;

      phi02d_own = aphi2v_own*rho02Sc_;
      phi02d_nei = aphi2v_nei*rho02Sc_;

    }


    void Foam::interTwoPhaseCentralFoam::ThermalConductivity()
    {

      alpha1_ = mu1_/Pr1_;
      alpha2_ = mu2_/Pr2_;

    }

/*
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

        ThermalConductivity();
        Tviscosity1 =
        (
          - fvm::laplacian(alpha1_*Cp1_, T_)
          - 0.5*fvc::laplacian(alpha1_, magSqr(U_))
        );

        Tviscosity2 =
        (
          - fvm::laplacian(alpha2_*Cp2_, T_)
          - 0.5*fvc::laplacian(alpha2_, magSqr(U_))
        );

    }


    void Foam::interTwoPhaseCentralFoam::devRhoReff()
    {

      ThermalConductivity();
      devRhoReff1_ = (-(alpha1_)*dev(twoSymm(fvc::grad(U_))));
      devRhoReff2_ = (-(alpha2_)*dev(twoSymm(fvc::grad(U_))));

      TSource1_ =
      fvc::div((linearInterpolate((-devRhoReff1_) & U_) & U_.mesh().Sf())());

      TSource2_ =
      fvc::div((linearInterpolate((-devRhoReff2_) & U_) & U_.mesh().Sf())());

    }
*/
    void Foam::interTwoPhaseCentralFoam::TSource()
    {

      rho1AmpU_ = magSqr(U_)*rho1_;
      rho2AmpU_ = magSqr(U_)*rho2_;

      TSource1_ = fvc::ddt(p_) - 0.5*fvc::ddt(rho1AmpU_);
      TSource2_ = fvc::ddt(p_) - 0.5*fvc::ddt(rho2AmpU_);
    }


    void Foam::interTwoPhaseCentralFoam::info()
    {

      E_ = volumeFraction1_*E1_ + volumeFraction2_*E2_;

      Info<< " max: E1 " << max(E1_).value()
          << " min: " << min(E1_).value()
          << " max: E2 " << max(E2_).value()
          << " min: " << min(E2_).value()
          << " max: E " << max(E_).value()
          << " min: " << min(E_).value()
          << nl << endl;

    }


    void Foam::interTwoPhaseCentralFoam::divU()
    {

      /*
              surfaceScalarField rbyAf = fvc::interpolate(rbyA_);

              divU();
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
