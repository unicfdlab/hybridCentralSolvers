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
#include "gaussGrad.H"

void Foam::vofTwoPhaseCentralFoam::UpdateCentralWeights()
{
    const auto& Sf = U_.mesh().Sf();
    const auto& magSf = U_.mesh().magSf();
    // phiv_own_ = fvc::interpolate(U_, own_, "reconstruct(U)") & Sf;
    // phiv_nei_ = fvc::interpolate(U_, nei_, "reconstruct(U)") & Sf;

    const auto &rho1 = mixture_model_.rho1();
    const auto &rho2 = mixture_model_.rho2();
    Info << "rho_phi"<< endl;
    surfaceScalarField rho_phi_own
    (
        (fvc::interpolate(rho1*U_, own_, "reconstruct(U)") & Sf)*vF1face_
        +
        (fvc::interpolate(rho2*U_, own_, "reconstruct(U)") & Sf)*vF2face_
    );
    surfaceScalarField rho_phi_nei
    (
        (fvc::interpolate(rho1*U_, nei_, "reconstruct(U)") & Sf)*vF1face_
        +
        (fvc::interpolate(rho2*U_, nei_, "reconstruct(U)") & Sf)*vF2face_
    );
    Info << "rho"<< endl;
    surfaceScalarField rho_own (vF1face_*rho1_own_ + vF2face_*rho2_own_);
    surfaceScalarField rho_nei (vF1face_*rho1_nei_ + vF2face_*rho2_nei_);
    phiv_own_ = rho_phi_own / rho_own;
    phiv_own_ = rho_phi_nei / rho_nei;
    Info << "CfSf"<< endl;
    CfSf_own_     = Cf_own_ * magSf;
    CfSf_own_.setOriented(true);
    CfSf_nei_     = Cf_nei_ * magSf;
    CfSf_nei_.setOriented(true);

    surfaceScalarField ap
    (
        max(max(phiv_own_ + CfSf_own_, phiv_nei_ + CfSf_nei_), v_zero_)
    );
    surfaceScalarField am
    (
        min(min(phiv_own_ - CfSf_own_, phiv_nei_ - CfSf_nei_), v_zero_)
    );

    alpha_own_   = ap/(ap - am);
    aSf_     = am*alpha_own_;
    alpha_nei_   = 1.0 - alpha_own_;
}

void Foam::vofTwoPhaseCentralFoam::UpdateCentralFields()
{
    const surfaceVectorField& Sf = U_.mesh().Sf();
    Info << "rAU" << endl;
    rAUf_own_ = alpha_own_*fvc::interpolate(oneByA_, own_, "reconstruct(rAU)");
    rAUf_nei_ = alpha_nei_*fvc::interpolate(oneByA_, nei_, "reconstruct(rAU)");
    Info << "phiHbyA" << endl;
    Info << "max/min HbyA:" << max(HbyA_) << "/" << min(HbyA_) << endl;
    phiHbyA_own_ =
        alpha_own_*((fvc::interpolate(HbyA_, own_, "reconstruct(U)")) & Sf)
        - aSf_;
    phiHbyA_nei_ =
        alpha_nei_*((fvc::interpolate(HbyA_, nei_, "reconstruct(U)")) & Sf)
        + aSf_;
    Info<<"done"<<endl;

    //This interpolation causes instability:
    // const auto &rho1 = mixture_model_.rho1();
    // const auto &rho2 = mixture_model_.rho2();
    // surfaceScalarField rho_phi_own =
    //     (fvc::interpolate(rho1*HbyA_, own_, "reconstruct(U)") & Sf)*vF1face_
    //     +
    //     (fvc::interpolate(rho2*HbyA_, own_, "reconstruct(U)") & Sf)*vF2face_;
    // surfaceScalarField rho_phi_nei =
    //     (fvc::interpolate(rho1*HbyA_, nei_, "reconstruct(U)") & Sf)*vF1face_
    //     +
    //     (fvc::interpolate(rho2*HbyA_, nei_, "reconstruct(U)") & Sf)*vF2face_;
    // surfaceScalarField rho_own =
    //     vF1face_*rho1_own_ + vF2face_*rho2_own_;
    // surfaceScalarField rho_nei =
    //     vF1face_*rho1_nei_ + vF2face_*rho2_nei_;
    // phiHbyA_own_ = rho_phi_own / rho_own - aSf_;
    // phiHbyA_own_ = rho_phi_nei / rho_nei + aSf_;
}

void Foam::vofTwoPhaseCentralFoam::CalculateMassFluxes()
{
    phi1_own_ = rho1_own_*aphiv_own_;
    phi1_nei_ = rho1_nei_*aphiv_nei_;

    phi2_own_ = rho2_own_*aphiv_own_;
    phi2_nei_ = rho2_nei_*aphiv_nei_;
}
//
//END-OF-FILE
//
