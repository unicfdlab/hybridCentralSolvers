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
#include "MULES.H"

//Equations of the model

void Foam::interTwoPhaseCentralFoam::solveRho
(
    volScalarField& rhoi,
    const surfaceScalarField& phii_own,
    const surfaceScalarField& phii_nei
)
{
    solve
    (
        fvm::ddt(rhoi)
        +
        fvc::div(phii_own)
        +
        fvc::div(phii_nei)
    );
}

void Foam::interTwoPhaseCentralFoam::solveRho1()
{
    solveRho(rho1_, phi1_own_, phi1_nei_);
    rho1_ = max(rho1_,rho1Min);
}

void Foam::interTwoPhaseCentralFoam::solveRho2()
{
    solveRho(rho2_, phi2_own_, phi2_nei_);
    rho2_ = max(rho2_,rho2Min);
}

void Foam::interTwoPhaseCentralFoam::alpha1Eqnsolve()
{
    //MULES version with antidiffusive flux

    const fvMesh& mesh = phi_.mesh();
    word alphaScheme("div(phi,volumeFraction1)");
    word alpharScheme("div(phirb,volumeFraction1)");

    interface_.correct();

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;

    // Set the time blending factor, 1 for Euler
    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    scalar cAlpha  =  interface_.cAlpha(); //or read from dictionary

    // Standard face-flux compression coefficient
    surfaceScalarField phic(cAlpha*mag(phi_/mesh.magSf()));

    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();

    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi_);

    volScalarField divU
    (
        fvc::div(phi_)
    );

    const dictionary& vf1Controls = mesh.solverDict(volumeFraction1_.name());
    label nAlphaCorr(vf1Controls.get<label>("nAlphaCorr"));

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                mesh.time().timeName(),
                mesh
            ),
            divU*0.0
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                mesh.time().timeName(),
                mesh
            ),
            (1.0 + K_)*divU*volumeFraction1_
        );


        surfaceScalarField phir(phic*interface_.nHatf());

        tmp<surfaceScalarField> talphaPhi1Un
        (
            fvc::flux
            (
                phiCN(),
                cnCoeff*volumeFraction1_,
                alphaScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, volumeFraction2_, alpharScheme),
                volumeFraction1_,
                alpharScheme
            )
        );

        surfaceScalarField alphaPhi10 = talphaPhi1Un;

        MULES::explicitSolve
        (
            geometricOneField(),
            volumeFraction1_,
            phiCN,
            alphaPhi10,
            Sp,
            Su,
            oneField(),
            zeroField()
        );

        //#include "alphasMinMax.H"
    }


     vF1face_ = fvc::interpolate
     (
         volumeFraction1_,
         "reconstruct(volumeFraction1)"
     );
     vF2face_ = 1.0 - vF1face_;
/*
     fvScalarMatrix alpha1Eqn
     (

         fvm::ddt(rho1_,volumeFraction1_)
         +
         fvm::div((phi1_own_ + phi1_nei_), volumeFraction1_)
     );

     alpha1Eqn.solve();
*/
     volumeFraction2_ = 1 - volumeFraction1_;

     Info<< "Phase-1 volume fraction = "
         << volumeFraction1_.weightedAverage(U_.mesh().Vsc()).value()
         << "  Min(" << volumeFraction1_.name() << ") = " << min(volumeFraction1_).value()
         << "  Max(" << volumeFraction1_.name() << ") = " << max(volumeFraction1_).value()
         << endl;
}


void Foam::interTwoPhaseCentralFoam::UEqn()
{
    Density();
    divDevRhoReff();

    surfaceScalarField phiU_own = vF1face_*phi1_own_ + vF2face_*phi2_own_;
    surfaceScalarField phiU_nei = vF1face_*phi1_nei_ + vF2face_*phi2_nei_;

    E_ = fvc::ddt(rho_) + fvc::div(phiU_own) + fvc::div(phiU_nei);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho_,U_) - fvm::Sp(E_,U_)
        +
        fvm::div(phiU_own,U_) + fvm::div(phiU_nei,U_)
        +
        divDevRhoReff_
    );

    rbyA_  = 1.0/UEqn.A();
    HbyA_ = UEqn.H()*rbyA_;
    HbyA_.boundaryFieldRef() == U_.boundaryField();
}

void Foam::interTwoPhaseCentralFoam::ReconstructVelocity()
{
    pressureGradient();
    U_ = HbyA_ - rbyA_*gradp_
         +  rbyA_*fvc::reconstruct(phib_/rAUf_);
    U_.correctBoundaryConditions();
}


void Foam::interTwoPhaseCentralFoam::TEqnsolve()
{
    fvScalarMatrix TEqn1
    (
        fvm::ddt(rho1_,T_)
        + fvm::div(phi1_own_,T_) + fvm::div(phi1_nei_,T_)
        - fvm::Sp(E1_,T_)
        + 1/Cp1_*TSource1_
    );

    fvScalarMatrix TEqn2
    (
        fvm::ddt(rho2_,T_)
        + fvm::div(phi2_own_,T_) + fvm::div(phi2_nei_,T_)
        - fvm::Sp(E2_,T_)
        + 1/Cp2_*TSource2_
    );

    fvScalarMatrix TEqn
    (
        T_,
        rho1_.dimensions()*T_.dimensions()/dimTime
    );

    combineMatrices
    (
        TEqn1,
        TEqn2,
        volumeFraction1_,
        volumeFraction2_,
        TEqn,
        true //do not copy: update origin matrices
    );

    TEqn.solve();
}


void Foam::interTwoPhaseCentralFoam::pEqnsolve()
{
    Wp_ = 1/(1 - (volumeFraction1_*psi1_ + volumeFraction2_*psi2_)*gh_);

//    p_rgh_ = (p_/Wp_) - rho0_*gh_;

    pEqn1_own_ =
    (
        fvc::div(phi01d_own_)
        + fvm::div(phi1d_own_,p_rgh_)
        - fvm::laplacian(Dp1_own_, p_rgh_)
    );

    pEqn1_nei_ =
    (
        fvc::div(phi01d_nei_)
        + fvm::div(phi1d_nei_,p_rgh_)
        - fvm::laplacian(Dp1_nei_, p_rgh_)
    );

    pEqn2_own_ =
    (
        fvc::div(phi02d_own_)
        + fvm::div(phi2d_own_,p_rgh_)
        - fvm::laplacian(Dp2_own_, p_rgh_)
    );

    pEqn2_nei_ =
    (
        fvc::div(phi02d_nei_)
        + fvm::div(phi2d_nei_,p_rgh_)
        - fvm::laplacian(Dp2_nei_, p_rgh_)
    );

    fvScalarMatrix pEqn1
    (
//        fvc::ddt(rho1_) + psi1_*correction(fvm::ddt(p_rgh_)) +
        fvm::ddt(psi1_,p_rgh_) +
        pEqn1_own_ + pEqn1_nei_
    );

    fvScalarMatrix pEqn2
    (
//        fvc::ddt(rho2_) + psi2_*correction(fvm::ddt(p_rgh_)) +
        fvm::ddt(psi2_,p_rgh_) +
        pEqn2_own_ + pEqn2_nei_
    );

    fvScalarMatrix pEqn
    (
        p_rgh_,
        p_rgh_.dimensions()*psi1_.dimensions()/dimTime
    );

    combineMatrices
    (
        pEqn1,
        pEqn2,
        volumeFraction1_,
        volumeFraction2_,
        pEqn,
        true //do not copy: update origin matrices
    );

    pEqn.solve();

    p_ = Wp_*(p_rgh_ + rho0_*gh_);

    const Foam::Time& runTime = U_.mesh().time();

    if(runTime.outputTime())
    {
      surfaceScalarField phi1
      ( "phi1", phi1_own_+phi1_nei_);

      surfaceScalarField phi2
      ( "phi2", phi2_own_+phi2_nei_);

      phi1.write();
      phi2.write();
    }
}

//
//END-OF-FILE
//
