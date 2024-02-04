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
#include "MULES.H"
#include "constrainPressure.H"

//Equations of the model

void Foam::vofTwoPhaseCentralFoam::solveRho
(
    volScalarField& rhoi,
    const surfaceScalarField& phii_own,
    const surfaceScalarField& phii_nei
)
{
    surfaceScalarField rhoiPhi = 
        phii_own + phii_nei;
    solve
    (
        fvm::ddt(rhoi)
        +
        fvc::div(rhoiPhi)
    );
}

void Foam::vofTwoPhaseCentralFoam::solveRho1()
{
    auto& rho1 = mixture_model_.rho1();
    const auto cdt = 0.0/rho1.mesh().time().deltaT();
    solveRho(rho1, phi1_own_, phi1_nei_);
    mixture_model_.limitRho1();
    E1_ = cdt * rho1;
}

void Foam::vofTwoPhaseCentralFoam::solveRho2()
{
    auto& rho2 = mixture_model_.rho2();
    const auto cdt = 0.0/rho2.mesh().time().deltaT();
    solveRho(rho2, phi2_own_, phi2_nei_);
    mixture_model_.limitRho2();
    E2_ = cdt * rho2;
}

void Foam::vofTwoPhaseCentralFoam::LiquidVolumeFractionSolve()
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
            (1.0 + Lambda_)*divU*volumeFraction1_
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

        // #include "alphasMinMax.H"
    }

    volumeFraction2_ = 1 - volumeFraction1_;

    dotVF1_ = fvc::ddt(volumeFraction1_);
    dotVF2_ = fvc::ddt(volumeFraction2_);

    vF1face_ = fvc::interpolate
    (
        volumeFraction1_,
        "reconstruct(volumeFraction1)"
    );
    vF2face_ = 1.0 - vF1face_;

    phiVF1_ = phi_*vF1face_;
    phiVF2_ = phi_ - phiVF1_;

    Info<< "Phase-1 volume fraction = "
        << volumeFraction1_.weightedAverage(U_.mesh().Vsc()).value()
        << "  Min(" << volumeFraction1_.name() << ") = " << min(volumeFraction1_).value()
        << "  Max(" << volumeFraction1_.name() << ") = " << max(volumeFraction1_).value()
        << endl;
}


void Foam::vofTwoPhaseCentralFoam::UEqn()
{
    const auto& rho = mixture_model_.rho();

    surfaceScalarField phiU_own = vF1face_*phi1_own_ +
        vF2face_*phi2_own_;
    surfaceScalarField phiU_nei = vF1face_*phi1_nei_ +
        vF2face_*phi2_nei_;
    phiU_own.rename("phiU_own");
    phiU_nei.rename("phiU_nei");

    E_ = fvc::ddt(rho) + fvc::div(phiU_own) + fvc::div(phiU_nei);
    fvVectorMatrix UEqn
    (
          fvm::ddt(rho,U_)
        + fvm::div(phiU_own,U_) + fvm::div(phiU_nei,U_)
        - fvm::Sp(E_,U_)
        + turbulence_->divDevRhoReff(rho,U_)
    );

    // const auto& rho1 = mixture_model_.rho1();
    // const auto& rho2 = mixture_model_.rho2();
    // fvVectorMatrix UEqn1
    // (
    //     fvm::ddt(rho1,U_)
    //     + fvm::div(phi1_own_,U_) + fvm::div(phi1_nei_,U_)
    //     - fvm::Sp(E1_,U_)
    // );
    // fvVectorMatrix UEqn2
    // (
    //     fvm::ddt(rho2,U_)
    //     + fvm::div(phi2_own_,U_) + fvm::div(phi2_nei_,U_)
    //     - fvm::Sp(E2_,U_)
    // );
    // fvVectorMatrix UEqn
    // (
    //     U_,
    //     UEqn1.dimensions()
    // );
    // combineMatrices
    // (
    //     UEqn1,
    //     UEqn2,
    //     volumeFraction1_,
    //     volumeFraction2_,
    //     UEqn,
    //     true //do not copy: update origin matrices
    // );
    // UEqn += turbulence_->divDevRhoReff(rho,U_);

    oneByA_  = 1.0/UEqn.A();
    HbyA_ = UEqn.H()*oneByA_;
    HbyA_.boundaryFieldRef() == U_.boundaryField();

    B_ = fvc::reconstruct
    (
        (
           interface_.surfaceTensionForce()
            - ghf_*fvc::snGrad(rho)
        ) * B_.mesh().magSf()
    );
}

void Foam::vofTwoPhaseCentralFoam::TEqnSolve()
{
    const auto& rho1 = mixture_model_.rho1();
    const auto& rho2 = mixture_model_.rho2();
    const auto& rho  = mixture_model_.rho();
    const auto& alpha1=mixture_model_.alpha1();
    const auto& alpha2=mixture_model_.alpha2();
    const auto& Cp1  = mixture_model_.Cp1();
    const auto& Cp2  = mixture_model_.Cp2();
    const auto& Prt  = mixture_model_.Prt();

    //yi = vFi*rhoi / rho
    
    volScalarField d1 = rho1 / rho;
    volScalarField d2 = rho2 / rho;

    volScalarField BdotU
    (
        B_ & U_
    );

    volScalarField Vd1
    (
        fvc::div //viscous dissipation for the first phase
        (
            (
                linearInterpolate
                (
                    (-rho1*turbulence_->devRhoReff()) & U_
                ) & U_.mesh().Sf()
            )()
        )
    );
    volScalarField Vd2
    (
        fvc::div //viscous dissipation for the second phase
        (
            (
                linearInterpolate
                (
                    (-rho2*turbulence_->devRhoReff()) & U_
                ) & U_.mesh().Sf()
            )()
        )
    );

    surfaceScalarField alphaEff=
        linearInterpolate(rho1*turbulence_->nuEff())/Prt
        +
        vF1face_*alpha1
        +
        vF2face_*alpha2;
    alphaEff.rename("alphaEff");

    fvScalarMatrix TEqn1
    (
        fvm::ddt(rho1,T_)
        + fvm::div(phi1_own_,T_) + fvm::div(phi1_nei_,T_)
        - fvm::Sp(E1_,T_)
        + 1.0/Cp1*TSource1_
        ==
        d1*BdotU/Cp1
        +
        Vd1/Cp1
    );

    fvScalarMatrix TEqn2
    (
        fvm::ddt(rho2,T_)
        + fvm::div(phi2_own_,T_) + fvm::div(phi2_nei_,T_)
        - fvm::Sp(E2_,T_)
        + 1.0/Cp2*TSource2_
        ==
        d2*BdotU/Cp2
        +
        Vd2/Cp2
    );

    fvScalarMatrix TEqn
    (
        T_,
        TEqn1.dimensions()
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

    solve
    (
        TEqn
        ==
        fvm::laplacian(alphaEff, T_)
    );
}


void Foam::vofTwoPhaseCentralFoam::pEqnSolve()
{
    const auto& rho  = mixture_model_.rho();
    const auto& rho0 = mixture_model_.rho0();
    const auto& rho1 = mixture_model_.rho1();
    const auto& rho2 = mixture_model_.rho2();
    const auto& psi1 = mixture_model_.psi1();
    const auto& psi2 = mixture_model_.psi2();

    phiHbyA_ = fvc::flux(HbyA_);
    surfaceScalarField rAUf("rAUf", linearInterpolate(oneByA_));

    phib_ = 
    (
        interface_.surfaceTensionForce()
        -
        ghf_*fvc::snGrad(rho)
    )* rAUf * rho.mesh().magSf();

    // blend KNP coeffs with linear interpolation
    rAUf_own_ *= kappa_; rAUf_own_ += onemkappa_*rAUf*alpha_own_;
    rAUf_nei_ *= kappa_; rAUf_nei_ += onemkappa_*rAUf*alpha_nei_;
    phiHbyA_own_ *= kappa_; phiHbyA_own_ += alpha_own_*onemkappa_*phiHbyA_;
    phiHbyA_nei_ *= kappa_; phiHbyA_nei_ += alpha_nei_*onemkappa_*phiHbyA_;

    phiHbyA_own_ += alpha_own_*onemkappa_*linearInterpolate(rho*oneByA_)*
        fvc::ddtCorr(U_, phi_);
    phiHbyA_nei_ += alpha_nei_*onemkappa_*linearInterpolate(rho*oneByA_)*
        fvc::ddtCorr(U_, phi_);
    phiHbyA_own_ += alpha_own_*phib_;
    phiHbyA_nei_ += alpha_nei_*phib_;

    // Update the pressure BCs to ensure flux consistency
    phiHbyA_      = phiHbyA_own_ + phiHbyA_nei_;
    constrainPressure(p_rgh_, U_, phiHbyA_, rAUf);

    // DDt(alpha1*rho1) =
    // ddt(alpha1*rho1) +
    // div(alpha1*rho1*U) -
    // rho1*(ddt(alpha1) + div(alpha1*U))
    // asymm fluxes:
    // ddt(alpha1*rho1) + fvc::div(alpha1*rho1*U) \approx
    // alpha1*(ddt(rho1) + fvc::div(rho1*U))
    // c = alpha1/rho1
    surfaceScalarField phiRho1 =
        aphiv_own_*rho1_own_ + aphiv_nei_*rho1_nei_;
    fvScalarMatrix cDDtAlpha1Rho1 = pos(volumeFraction1_)*
    (
        (volumeFraction1_/rho1)*
        (
            fvc::ddt(rho1) + fvc::div(phiRho1) + psi1*correction(fvm::ddt(p_rgh_))
        )
        - dotVF1_ - fvc::div(phiVF1_)
    );

    surfaceScalarField phiRho2 =
        aphiv_own_*rho2_own_ + aphiv_nei_*rho2_nei_;
    fvScalarMatrix cDDtAlpha2Rho2 = pos(volumeFraction2_)*
    (
        (volumeFraction2_/rho2)*
        (
            fvc::ddt(rho2) + fvc::div(phiRho2) + psi2*correction(fvm::ddt(p_rgh_))
        )
        - dotVF2_ - fvc::div(phiVF2_)
    );

    fvScalarMatrix elliptic_p_own = fvc::div(phiHbyA_own_) - fvm::laplacian(rAUf_own_,p_rgh_);
    fvScalarMatrix elliptic_p_nei = fvc::div(phiHbyA_nei_) - fvm::laplacian(rAUf_nei_,p_rgh_);

    solve
    (
        cDDtAlpha1Rho1 +
        cDDtAlpha2Rho2 +
        elliptic_p_own +
        elliptic_p_nei
    );
    
    aphiv_own_ = phiHbyA_own_ + elliptic_p_own.flux();
    aphiv_nei_ = phiHbyA_nei_ + elliptic_p_nei.flux();
    phi_ = aphiv_own_ + aphiv_nei_;

    U_ = HbyA_
    //    + oneByA_*
    //    (
    //         fvc::reconstruct
    //         (
    //             phib_/
    //             (rAUf_own_ + rAUf_nei_)
    //         )
    //         -
    //         fvc::grad(p_rgh_)
    //    );
       + oneByA_*
       fvc::reconstruct
       (
           (phib_ + elliptic_p_own.flux() + elliptic_p_nei.flux())/
           (rAUf_own_ + rAUf_nei_)
       );
    U_.correctBoundaryConditions();
    //fvOptions.correct(U_);


    Wp_ = 1.0/(1.0 - (volumeFraction1_*psi1 + volumeFraction2_*psi2)*gh_);

    p_ = Wp_*(p_rgh_ + rho0*gh_);
    p_.correctBoundaryConditions();
}

//
//END-OF-FILE
//
