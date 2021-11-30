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
}

void Foam::interTwoPhaseCentralFoam::solveRho2()
{
    solveRho(rho2_, phi2_own_, phi2_nei_);
}

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
            - fvm::Sp(E1_,T_)
//          + Tviscosity1
            + 1/Cp1_*TSource1_
        )
        + volumeFraction2_*
        (
            fvm::ddt(rho2_,T_)
            + fvm::div(phi2_own_,T_) + fvm::div(phi2_nei_,T_)
            - fvm::Sp(E2_,T_)
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

//
//END-OF-FILE
//

