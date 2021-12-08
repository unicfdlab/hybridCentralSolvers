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

    interface_
    (
        volumeFraction1_,
        U_,
        *this
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

    Cf_own_
    (
        "Cf_own",
        0.0*aSf1_/mesh.magSf()
    ),

    Cf_nei_
    (
        "Cf_nei",
        Cf_own_
    ),

    mu_
    (
        "mu",
        volumeFraction1_*mu1_
        +
        volumeFraction2_*mu2_
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

    divDevRhoReff_
    (
        - fvm::laplacian(mu1_, U_)
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
    rho_.oldTime();
    U_.oldTime();
    T_.oldTime();
    p_.oldTime();
    psi1_.oldTime();
    psi2_.oldTime();
    Q_.oldTime();
}

void Foam::interTwoPhaseCentralFoam::massError1()
{
    E1_ =
    (
        fvc::ddt(rho1_) + fvc::div(phi1_own_ + phi1_nei_)
    );

    Info<< " max: rho1 " << max(rho1_).value()
        << " min: " << min(rho1_).value()
        << endl;
}


void Foam::interTwoPhaseCentralFoam::massError2()
{
    E2_ =
    (
        fvc::ddt(rho2_) + fvc::div(phi2_own_ + phi2_nei_)
    );

    Info<< " max: rho2 " << max(rho2_).value()
        << " min: " << min(rho2_).value()
        << endl;
}


void Foam::interTwoPhaseCentralFoam::TSource()
{
    Q_ = 0.5*magSqr(U_);

    TSource1_ =
    (
      fvc::ddt(rho1_,Q_)
      + fvc::div(phi1_own_,Q_) + fvc::div(phi1_nei_,Q_)
      - fvc::ddt(p_)
      - fvc::Sp(E1_,Q_)
    );

    TSource2_ =
    (
      fvc::ddt(rho2_,Q_)
      + fvc::div(phi2_own_,Q_) + fvc::div(phi2_nei_,Q_)
      - fvc::ddt(p_)
      - fvc::Sp(E2_,Q_)
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

    updateK();

    // Thermal conductivity
    alpha1_ = mu1_/Pr1_;
    alpha2_ = mu2_/Pr2_;

    UpdateCentralWeightsIndividual();

//      divDevRhoReff();

    UEqn();

    UpdateCentralFieldsIndividual();

    pressureGradient();

//      devRhoReff();

//      Tviscosity1 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_));
//      Tviscosity2 = (- 0*fvm::laplacian(alpha1_*Cp1_, T_));

    TSource1_ = 0*fvc::ddt(p_);
    TSource2_ = 0*fvc::ddt(p_);
    TSource_ = volumeFraction1_*1/Cp1_*TSource1_;
}

//* * * * * * * * * * * * * * * * * Flux` Functions * * * * * * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::Flux()
{
    phi1_own_ = pEqn1_own_.flux() + phi01d_own_;
    phi1_nei_ = pEqn1_nei_.flux() + phi01d_nei_;

    phi2_own_ = pEqn2_own_.flux() + phi02d_own_;
    phi2_nei_ = pEqn2_nei_.flux() + phi02d_nei_;
}

void Foam::interTwoPhaseCentralFoam::updateKappa()
{
    const fvMesh& mesh = U_.mesh();
    surfaceScalarField phiv_own = fvc::interpolate
        (
            U_,
            own_,
            "reconstruct(U)"
        ) & mesh.Sf();
    surfaceScalarField phiv_nei = fvc::interpolate
        (
            U_,
            nei_,
            "reconstruct(U)"
        ) & mesh.Sf();

    surfaceScalarField CfSf = max(Cf_own_, Cf_nei_) * mesh.magSf();
    CfSf.setOriented(true);

    surfaceScalarField amaxSf =
        max
        (
            max
            (
                mag(phiv_own + CfSf),
                mag(phiv_own - CfSf)
            ),
            max
            (
                mag(phiv_nei + CfSf),
                mag(phiv_nei - CfSf)
            )
        );
    amaxSf.setOriented(true);

    surfaceScalarField amaxSfbyDelta
    (
        mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
    );

    surfaceScalarField FaceAcCo =
        (
            amaxSfbyDelta/mesh.magSf() * mesh.time().deltaT()
        );

    surfaceScalarField Maf = linearInterpolate(mag(U_)/C_);

    kappa_ =
        min
        (
            Maf/FaceAcCo,
            scalar(1.0)
        );

    Info<< "max/min kappa: " << max(kappa_).value()
        << "/" << min(kappa_).value()
        << endl;

    writeMaxMinKappa (kappa_);

    //kappaBlend(kappa_, phi1_own_, phi1_nei_);
    //kappaBlend(kappa_, phi2_own_, phi2_nei_);
}

void Foam::interTwoPhaseCentralFoam::writeMaxMinKappa
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


void Foam::interTwoPhaseCentralFoam::volumeFlux()
{
    surfaceScalarField rho1_own = fvc::interpolate
    (
        rho1_,
        own_,
        "reconstruct(rho1)"
    );
    surfaceScalarField rho1_nei = fvc::interpolate
    (
        rho1_,
        nei_,
        "reconstruct(rho1)"
    );

    phi_ =
    (
        (phi1_own_ + phi1_nei_)
        /(alpha1_own_*rho1_own + alpha1_nei_*rho1_nei)
    );
}

void Foam::interTwoPhaseCentralFoam::cellCourantNo()
{
    // autoPtr<volScalarField> localCoPtr;
    // if (mesh.nInternalFaces())
    // {
    //     localCoPtr.reset
    //     (
    //         new volScalarField
    //         (
    //             IOobject
    //             (
    //                 "Co",
    //                 mesh.time().timeName(),
    //                 mesh
    //             ),
    //             mesh,
    //             dimensionedScalar("localCo", dimless,  VSMALL),
    //             zeroGradientFvPatchScalarField::typeName
    //         )
    //     );

    //     surfaceScalarField rhoPhi =

    //     localCoPtr().primitiveFieldRef() =
    //     (
    //         fvc::surfaceSum(mag(phi))().primitiveField()
    //         /
    //         rho.primitiveField()
    //         /
    //         mesh.V().field()
    //     ) * 0.5 * runTime.deltaT().value();

    //     CoNum = gMax(localCoPtr().internalField());
    // }

    // Info<< "Courant Number max: " <<  CoNum << endl;
}

void Foam::interTwoPhaseCentralFoam::combineMatrices
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

void Foam::interTwoPhaseCentralFoam::updateK()
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

    Cf_own_ = fvc::interpolate(C_, own_, "reconstruct(psi)");
    Cf_nei_ = fvc::interpolate(C_, nei_, "reconstruct(psi)");
}

void Foam::interTwoPhaseCentralFoam::UpdateCentralWeights
(
    const volScalarField& rhoi,
    const volVectorField& U,
    const volScalarField& Ci,
    surfaceScalarField& alpha_own,
    surfaceScalarField& alpha_nei,
    surfaceScalarField& aSf
)
{
    const auto& Sf = U.mesh().Sf();
    const auto&mSf = U.mesh().magSf();

    surfaceScalarField rhoi_own =
        fvc::interpolate(rhoi, own_, "reconstruct(" + rhoi.name() + ")");
    surfaceScalarField rhoi_nei =
        fvc::interpolate(rhoi, nei_, "reconstruct(" + rhoi.name() + ")");

    surfaceScalarField phiv_own =
        (fvc::interpolate(rhoi*U, own_, "reconstruct(U)") & Sf)/rhoi_own;
    surfaceScalarField phiv_nei =
        (fvc::interpolate(rhoi*U, nei_, "reconstruct(U)") & Sf)/rhoi_nei;

    surfaceScalarField Ci_own = fvc::interpolate(Ci, own_, "reconstruct(psi)");
    surfaceScalarField Ci_nei = fvc::interpolate(Ci, nei_, "reconstruct(psi)");

    surfaceScalarField CiSf_own = Ci_own*mSf;
    CiSf_own.setOriented(true);
    surfaceScalarField CiSf_nei = Ci_nei*mSf;
    CiSf_nei.setOriented(true);

    surfaceScalarField api =
        max(max(phiv_own + CiSf_own, phiv_nei + CiSf_nei), v_zero);
    surfaceScalarField ami =
        min(min(phiv_own - CiSf_own, phiv_nei - CiSf_nei), v_zero);

    alpha_own = api/(api - ami);
    aSf = ami*alpha_own;
    alpha_nei = 1.0 - alpha_own;
}

void Foam::interTwoPhaseCentralFoam::UpdateCentralMassFluxes
(
    const volScalarField& rhoi,
    const dimensionedScalar& rho0i,
    const volScalarField& psii,
    const volVectorField& HbyA,
    const volScalarField& rbyA,
    const surfaceScalarField& alpha_own,
    const surfaceScalarField& alpha_nei,
    const surfaceScalarField& aSf,
    surfaceScalarField& phidi_own,
    surfaceScalarField& phidi_nei,
    surfaceScalarField& phi0i_own,
    surfaceScalarField& phi0i_nei,
    surfaceScalarField& Dpi_own,
    surfaceScalarField& Dpi_nei
)
{
    surfaceScalarField rhoi_own =
        fvc::interpolate(rhoi, own_, "reconstruct(" + rhoi.name() + ")");
    surfaceScalarField rhoi_nei =
        fvc::interpolate(rhoi, nei_, "reconstruct(" + rhoi.name() + ")");

    surfaceScalarField psii_own =
        fvc::interpolate(psii, own_, "reconstruct(" + psii.name() + ")");
    surfaceScalarField psii_nei =
        fvc::interpolate(psii, nei_, "reconstruct(" + psii.name() + ")");

    surfaceVectorField rhoUi_own =
        fvc::interpolate(rhoi*HbyA, own_, "reconstruct(U)");
    surfaceVectorField rhoUi_nei =
        fvc::interpolate(rhoi*HbyA, nei_, "reconstruct(U)");

    surfaceScalarField phiv_own = (rhoUi_own & U_.mesh().Sf())/rhoi_own;
    surfaceScalarField phiv_nei = (rhoUi_nei & U_.mesh().Sf())/rhoi_nei;

    surfaceScalarField aphiv_own = alpha_own*phiv_own - aSf;
    surfaceScalarField aphiv_nei = alpha_nei*phiv_nei + aSf;

    phidi_own = aphiv_own*psii_own;
    phidi_nei = aphiv_nei*psii_nei;

    phi0i_own = aphiv_own*rho0i;
    phi0i_nei = aphiv_nei*rho0i;

    Dpi_own = alpha_own*fvc::interpolate(rhoi*rbyA, own_, "reconstruct(Dp)");
    Dpi_nei = alpha_nei*fvc::interpolate(rhoi*rbyA, nei_, "reconstruct(Dp)");
}

void Foam::interTwoPhaseCentralFoam::kappaBlend
(
    const surfaceScalarField& kappa,
    surfaceScalarField& phi_own,
    surfaceScalarField& phi_nei
)
{
    phi_own += (1.0 - kappa) * phi_nei;
    phi_nei = kappa * phi_nei;
}

//* * * * * * * * * * * * Kurganov's coefficients Mixture * * * * * * * * * *//
void Foam::interTwoPhaseCentralFoam::UpdateCentralWeights()
{
    //Use uniform 1 density field to avoid density weighting
    volScalarField unitRho
    (
        "rho1",
        rho1_ / rho1_
    );
    speedOfSound();
    UpdateCentralWeights
    (
        unitRho,
        U_,
        C_,
        alpha_own_,
        alpha_nei_,
        aSf_
    );
}

void Foam::interTwoPhaseCentralFoam::UpdateCentralFields()
{
    // Update Phase 1
    UpdateCentralMassFluxes
    (
        rho1_,
        rho01_,
        psi1_,
        HbyA_,
        rbyA_,
        alpha_own_,
        alpha_nei_,
        aSf_,
        phi1d_own_,
        phi1d_nei_,
        phi01d_own_,
        phi01d_nei_,
        Dp1_own_,
        Dp1_nei_
    );

    // Update Phase 2
    UpdateCentralMassFluxes
    (
        rho2_,
        rho02_,
        psi2_,
        HbyA_,
        rbyA_,
        alpha_own_,
        alpha_nei_,
        aSf_,
        phi2d_own_,
        phi2d_nei_,
        phi02d_own_,
        phi02d_nei_,
        Dp2_own_,
        Dp2_nei_
    );
}

//* * * * * * * * * Kurganov's coefficients Individual Phases * * * * * * * *//

void Foam::interTwoPhaseCentralFoam::UpdateCentralWeightsIndividual()
{
//************************************ Phase One *****************************//

    volScalarField C1 = sqrt(gamma1_*R1_*T_);
    UpdateCentralWeights
    (
        rho1_,
        U_,
        C1,
        alpha1_own_,
        alpha1_nei_,
        aSf1_
    );

//************************************ Phase Two *****************************//

    volScalarField C2 = sqrt(gamma2_*R2_*T_);
    UpdateCentralWeights
    (
        rho2_,
        U_,
        C2,
        alpha2_own_,
        alpha2_nei_,
        aSf2_
    );

    Info<< " C1 " << max(C1).value()
        << " C2 " << max(C2).value()
        << endl;
}


void Foam::interTwoPhaseCentralFoam::UpdateCentralFieldsIndividual()
{
    // Update Phase 1
    UpdateCentralMassFluxes
    (
        rho1_,
        rho01_,
        psi1_,
        HbyA_,
        rbyA_,
        alpha1_own_,
        alpha1_nei_,
        aSf1_,
        phi1d_own_,
        phi1d_nei_,
        phi01d_own_,
        phi01d_nei_,
        Dp1_own_,
        Dp1_nei_
    );

    // Update Phase 2
    UpdateCentralMassFluxes
    (
        rho2_,
        rho02_,
        psi2_,
        HbyA_,
        rbyA_,
        alpha2_own_,
        alpha2_nei_,
        aSf2_,
        phi2d_own_,
        phi2d_nei_,
        phi02d_own_,
        phi02d_nei_,
        Dp2_own_,
        Dp2_nei_
    );
}
