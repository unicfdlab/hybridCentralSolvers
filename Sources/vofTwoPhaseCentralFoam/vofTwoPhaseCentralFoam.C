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

    R1_
    (
        dimensioned<scalar>("R1", *this)
    ),

    R2_
    (
        dimensioned<scalar>("R2", *this)
    ),

    molM1_
    (
        dimensioned<scalar>("molM1", *this)
    ),

    molM2_
    (
        dimensioned<scalar>("molM2", *this)
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
        1.0 - volumeFraction1_
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
        dimensioned<scalar>("rho01", *this)
    ),

    rho02_
    (
        dimensioned<scalar>("rho02", *this)
    ),

    Cp1_
    (
        dimensioned<scalar>("Cp1", *this)
    ),

    Cp2_
    (
        dimensioned<scalar>("Cp2", *this)
    ),

    mu1_
    (
        dimensioned<scalar>("mu1", *this)
    ),

    mu2_
    (
        dimensioned<scalar>("mu2", *this)
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
        dimensioned<scalar>("gamma1", Cp1_/(Cp1_ - (Foam::constant::physicoChemical::R/molM1_)))
    ),

    gamma2_
    (
        dimensioned<scalar>("gamma2", Cp2_/(Cp2_ - (Foam::constant::physicoChemical::R/molM2_)))
    ),

    C_
    (
        "C",
        sqrt(gamma1_*R1_*T_)
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

//     phib_
//     (
//         "phib",
//         0.0*(linearInterpolate(HbyA_) & mesh.Sf())
//     ),

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
        "phi1_own_",
        phi_*rho01_*0.0
    ),

    phi1_nei_
    (
        "phi1_nei_",
        phi_*rho01_*0.0
    ),

    phi2_own_
    (
        "phi2_own_",
        phi_*rho02_*0.0
    ),

    phi2_nei_
    (
        "phi2_nei_",
        phi_*rho02_*0.0
    ),

    rho1_own_
    (
        "rho1_own_",
        fvc::interpolate(rho1_, own_, "reconstruct(rho)")
    ),

    rho1_nei_
    (
        "rho1_nei_",
        fvc::interpolate(rho1_, nei_, "reconstruct(rho)")
    ),

    rho2_own_
    (
        "rho2_own_",
        fvc::interpolate(rho2_, own_, "reconstruct(rho)")
    ),

    rho2_nei_
    (
        "rho2_nei_",
        fvc::interpolate(rho2_, nei_, "reconstruct(rho)")
    ),

    alpha_own_
    (
        "alpha_own_ ",
        own_
    ),

    alpha_nei_
    (
        "alpha_nei_",
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

    mu_
    (
        "mu",
        volumeFraction1_*mu1_
        +
        volumeFraction2_*mu2_
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

    divDevRhoReff_
    (
        - fvm::laplacian(mu1_, U_)
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
        fvc::ddt(rho1_)
    ),

    E2_
    (
        fvc::ddt(rho2_)
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

    rho1Min_
    (
        dimensioned<scalar>("rho1Min", *this)
    ),

    rho2Min_
    (
        dimensioned<scalar>("rho2Min", *this)
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
        1.0/(1.0 - (volumeFraction1_*psi1_ + volumeFraction2_*psi2_)*gh_)
    ),

    rho0_
    (
        "rho0",
        volumeFraction1_*rho01_ + volumeFraction2_*rho02_
    ),

    B_
    (
        "B",
        0.0*rho0_*U_/U_.mesh().time().deltaT()
    ),

    phib_
    (
        "phib",
        0.0*phi_
    )

{
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
    rho1_.oldTime();
    rho2_.oldTime();
    rho_.oldTime();
    U_.oldTime();
    T_.oldTime();
    p_.oldTime();
    psi1_.oldTime();
    psi2_.oldTime();
    Q_.oldTime();
    p_rgh_.oldTime();
    phi_.oldTime();
}

void Foam::vofTwoPhaseCentralFoam::CharacteristicCourant()
{
    const auto &magSf = phi_.mesh().magSf();
    const auto &deltaCoeffs = phi_.mesh().deltaCoeffs();
    const auto &deltaT = phi_.mesh().time().deltaT();

    surfaceScalarField CfSf = max(CfSf_own_, CfSf_nei_);
    amaxSf_ = max(mag(phi_ + CfSf), mag(phi_ - CfSf));
    // surfaceScalarField uPlusC_pos =
    //     max(max(phi_own_ + CfSf_own_,phi_nei_ + CfSf_nei_),v_zero_)/magSf;
    // surfaceScalarField uPlusC_neg =
    //     min(min(phi_own_ - CfSf_own_,phi_nei_ - CfSf_nei_),v_zero_)/magSf;
    // surfaceScalarField uPlusC_max = max(uPlusC_pos,-uPlusC_neg);
    
    surfaceScalarField CCof = amaxSf_ * deltaCoeffs * deltaT / magSf;
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
    E1_ =
    (
        fvc::ddt(rho1_) + fvc::div(phi1_own_ + phi1_nei_)
    );

    Info<< " max: rho1 " << max(rho1_).value()
        << " min: " << min(rho1_).value()
        << endl;
}


void Foam::vofTwoPhaseCentralFoam::massError2()
{
    E2_ =
    (
        fvc::ddt(rho2_) + fvc::div(phi2_own_ + phi2_nei_)
    );

    Info<< " max: rho2 " << max(rho2_).value()
        << " min: " << min(rho2_).value()
        << endl;
}


void Foam::vofTwoPhaseCentralFoam::TSource()
{
    dpdt_ = (p_ - p_.oldTime()) / p_.mesh().time().deltaT();
    Q_ = 0.5*magSqr(U_);
    Q_.rename("Q");

    TSource1_ =
    (
        fvc::ddt(rho1_,Q_)
      + fvc::div(phi1_own_,Q_) + fvc::div(phi1_nei_,Q_)
      - dpdt_
      - fvc::Sp(E1_,Q_)
    );

    TSource2_ =
    (
        fvc::ddt(rho2_,Q_)
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

    Density();

    updateLambda();

    speedOfSound();

    interpolateDensities();

    CharacteristicCourant();

    // Thermal conductivity
    alpha1_ = mu1_/Pr1_;
    alpha2_ = mu2_/Pr2_;

    UpdateCentralWeights();
    UpdateCentralFields();

    UEqn();

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

    if (kappaIsOne)
    {
        kappa_.primitiveFieldRef() = 1.0;
        kappa_.boundaryFieldRef() = 1.0;
    }
    else
    {
        const fvMesh& mesh = U_.mesh();

        surfaceScalarField CfSf = max(CfSf_own_, CfSf_nei_);
        CfSf.setOriented(true);

        surfaceScalarField amaxSfbyDelta
        (
            mesh.surfaceInterpolation::deltaCoeffs()*amaxSf_
        );

        surfaceScalarField FaceAcCo =
        (
            amaxSfbyDelta/mesh.magSf() * mesh.time().deltaT()
        );


        surfaceScalarField Maf = mag(phi_) / CfSf;

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

//* * * * * * * * * * * * * * * Update Dencities * * * * * * * * * * * * * * *//

void Foam::vofTwoPhaseCentralFoam::DensityThermo()
{
    Info<< "Calculating phase densities"
        << endl;
    rho1_ = psi1_*p_ + rho01_;
    rho2_ = psi2_*p_ + rho02_;
    rho1_ = max(rho1_,rho1Min_);
    rho2_ = max(rho2_,rho2Min_);
    rho1_.correctBoundaryConditions();
    rho2_.correctBoundaryConditions();
}


void Foam::vofTwoPhaseCentralFoam::Density()
{
    Info<< "Calculating mixture density"
        << endl;
    rho_ = volumeFraction1_*rho1_ + volumeFraction2_*rho2_;
}

//* * * * * * * * * * * * * Update Dependable Variables * * * * * * * * * * *//

void Foam::vofTwoPhaseCentralFoam::updateLambda()
{
    Info<< "Calculating interface compressibility"
        << endl;
    volScalarField C1sqr = gamma1_*R1_*T_;
    volScalarField C2sqr = gamma2_*R2_*T_;

    volScalarField Z1 = rho1_*C1sqr;
    volScalarField Z2 = rho2_*C2sqr;

    Lambda_ =
    (
        (volumeFraction1_*volumeFraction2_*(Z2 - Z1))
        /(Z1*volumeFraction2_ + Z2*volumeFraction1_)
    );
}


void Foam::vofTwoPhaseCentralFoam::Compressibility()
{
    Info<< "Calculating phase compressibilities"
        << endl;
    psi1_ = 1.0/(R1_*T_);
    psi2_ = 1.0/(R2_*T_);
}


void Foam::vofTwoPhaseCentralFoam::speedOfSound()
{
    // volScalarField rbypsiM =
    //     1/(volumeFraction1_*psi1_ + volumeFraction2_*psi2_);
    // volScalarField psiM = 1/rbypsiM;

    volScalarField psiM = volumeFraction1_*psi1_ + volumeFraction2_*psi2_;

    volScalarField y1 = volumeFraction1_*(rho1_/rho_);

    volScalarField y2 = volumeFraction2_*(rho2_/rho_);

    volScalarField CpM = y1*Cp1_ + y2*Cp2_;

    volScalarField CvM = y1*(Cp1_/gamma1_) + y2*(Cp2_/gamma2_);

    volScalarField gammaM = CpM/CvM;

    C_ = sqrt(gammaM/psiM);

    Info<< "Calculating speed of sound interpolations"
        << endl;
    Cf_own_ = fvc::interpolate(C_, own_, "reconstruct(psi)");
    Cf_nei_ = fvc::interpolate(C_, nei_, "reconstruct(psi)");
}

void Foam::vofTwoPhaseCentralFoam::interpolateDensities()
{
    Info<< "Calculating phase densities interpolations"
        << endl;
    rho1_own_ =
        fvc::interpolate(rho1_, own_, "reconstruct(" + rho1_.name() + ")");
    rho1_nei_ =
        fvc::interpolate(rho1_, nei_, "reconstruct(" + rho1_.name() + ")");

    rho2_own_ =
        fvc::interpolate(rho2_, own_, "reconstruct(" + rho2_.name() + ")");
    rho2_nei_ =
        fvc::interpolate(rho2_, nei_, "reconstruct(" + rho2_.name() + ")");
}

void Foam::vofTwoPhaseCentralFoam::UpdateCentralWeights
(
    const volScalarField& rhoi,
    const surfaceScalarField& volumeFractioni,
    const volVectorField& U,
    const volScalarField& Ci,
    surfaceScalarField& rhoi_own,
    surfaceScalarField& rhoi_nei,
    surfaceScalarField& phiv_own,
    surfaceScalarField& phiv_nei,
    surfaceScalarField& alpha_own,
    surfaceScalarField& alpha_nei,
    surfaceScalarField& aSf
)
{
    const auto& Sf = U.mesh().Sf();
    const auto&mSf = U.mesh().magSf();

    rhoi_own =
        fvc::interpolate(rhoi, own_, "reconstruct(" + rhoi.name() + ")");
    rhoi_nei =
        fvc::interpolate(rhoi, nei_, "reconstruct(" + rhoi.name() + ")");

    phiv_own =
        volumeFractioni*(fvc::interpolate(rhoi*U, own_, "reconstruct(U)") & Sf)
        /rhoi_own;
    phiv_nei =
        volumeFractioni*(fvc::interpolate(rhoi*U, nei_, "reconstruct(U)") & Sf)
        /rhoi_nei;

    surfaceScalarField Ci_own = fvc::interpolate(Ci, own_, "reconstruct(psi)");
    surfaceScalarField Ci_nei = fvc::interpolate(Ci, nei_, "reconstruct(psi)");

    surfaceScalarField CiSf_own = Ci_own*mSf;
    CiSf_own.setOriented(true);
    surfaceScalarField CiSf_nei = Ci_nei*mSf;
    CiSf_nei.setOriented(true);

    surfaceScalarField api =
        max(max(phiv_own + CiSf_own, phiv_nei + CiSf_nei), v_zero_);
    surfaceScalarField ami =
        min(min(phiv_own - CiSf_own, phiv_nei - CiSf_nei), v_zero_);

    alpha_own = api/(api - ami);
    aSf = ami*alpha_own;
    alpha_nei = 1.0 - alpha_own;
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

//* * * * * * * * * Kurganov's coefficients Individual Phases * * * * * * * *//

void Foam::vofTwoPhaseCentralFoam::UpdateCentralWeights()
{
    Info<< "Calculating central weights"
        << endl;
    const auto& Sf = U_.mesh().Sf();
    const auto& magSf = U_.mesh().magSf();
    // take fluxes from the previous iteration (or the time step)
    phiv_own_ = fvc::interpolate(U_, own_, "reconstruct(U)") & Sf;
    phiv_nei_ = fvc::interpolate(U_, nei_, "reconstruct(U)") & Sf;

    CfSf_own_     = Cf_own_ * magSf;
    CfSf_own_.setOriented(true);
    CfSf_nei_     = Cf_nei_ * magSf;
    CfSf_nei_.setOriented(true);

    surfaceScalarField ap =
        max(max(phiv_own_ + CfSf_own_, phiv_nei_ + CfSf_nei_), v_zero_); //??? phiv_own ???
    surfaceScalarField am =
        min(min(phiv_own_ - CfSf_own_, phiv_nei_ - CfSf_nei_), v_zero_); //??? phiv_nei ???

    alpha_own_   = ap/(ap - am);
    aSf_     = am*alpha_own_;
    alpha_nei_   = 1.0 - alpha_own_;
}

void Foam::vofTwoPhaseCentralFoam::UpdateCentralFields()
{
    Info<< "Calculating central fields"
        << endl;
    const auto& Sf = U_.mesh().Sf();
	rAUf_own_ = alpha_own_*fvc::interpolate(oneByA_, own_, "reconstruct(rAU)");
	rAUf_nei_ = alpha_nei_*fvc::interpolate(oneByA_, nei_, "reconstruct(rAU)");

	phiHbyA_own_ =
        alpha_own_*((fvc::interpolate(HbyA_, own_, "reconstruct(U)")) & Sf)
        - aSf_;
	phiHbyA_nei_ =
        alpha_nei_*((fvc::interpolate(HbyA_, nei_, "reconstruct(U)")) & Sf)
        + aSf_;

    // GravityCorrection();

    // //add contribution from body forces
    // {
    //     surfaceScalarField ghSf = ghf_ * U_.mesh().magSf();

    //     rAUf_ = linearInterpolate(rbyA_);
    //     phib_ = -ghSf*fvc::snGrad(rho_)*rAUf_;

    //     {
    //         forAll(phib_.boundaryField(), patchi)
    //         {
    //             if (!phib_.boundaryField()[patchi].coupled())
    //             {
    //                 phib_.boundaryFieldRef()[patchi] = 0.0;
    //             }
    //         }
    //     }

    //     surfaceScalarField rho1f = linearInterpolate(rho1_);
    //     surfaceScalarField rho2f = linearInterpolate(rho2_);

    //     //surfaceScalarField phiCorr = linearInterpolate(rho_*rbyA_)*fvc::ddtCorr(U_, phi_);
    //     phi_ = linearInterpolate(HbyA_) & HbyA_.mesh().Sf();
    //     //phi_ += phiCorr;
    //     surfaceScalarField phikappa = phi_*onemkappa_;

    //     surfaceScalarField kapparAuf = onemkappa_*rAUf_;

    //     phi01d_own_ *= kappa_; phi01d_own_ += rho01_*phikappa;
    //     phi02d_own_ *= kappa_; phi02d_own_ += rho02_*phikappa;
    //     phi1d_own_  *= kappa_; phi1d_own_  += linearInterpolate(psi1_)*phikappa;
    //     phi2d_own_  *= kappa_; phi2d_own_  += linearInterpolate(psi2_)*phikappa;
    //     Dp1_own_    *= kappa_; Dp1_own_    += rho1f*kapparAuf;
    //     Dp2_own_    *= kappa_; Dp2_own_    += rho2f*kapparAuf;

    //     phi01d_nei_ *= kappa_;
    //     phi02d_nei_ *= kappa_;
    //     phi1d_nei_  *= kappa_;
    //     phi2d_nei_  *= kappa_;
    //     Dp1_nei_    *= kappa_;
    //     Dp2_nei_    *= kappa_;

    //     phi01d_own_ += rho1f*phib_;
    //     phi02d_own_ += rho2f*phib_;
    // }
}

void Foam::vofTwoPhaseCentralFoam::CalculateMassFluxes()
{
    phi1_own_ = rho1_own_*aphiv_own_;
    phi1_nei_ = rho1_nei_*aphiv_nei_;

    phi2_own_ = rho2_own_*aphiv_own_;
    phi2_nei_ = rho2_nei_*aphiv_nei_;
}

void Foam::vofTwoPhaseCentralFoam::GravityCorrection()
{
    // rho0_ = volumeFraction1_*rho01_ + volumeFraction2_*rho02_;

    // surfaceScalarField psi1_own =
    //     fvc::interpolate(psi1_, own_, "reconstruct(psi1)");

    // surfaceScalarField psi1_nei =
    //     fvc::interpolate(psi1_, nei_, "reconstruct(psi1)");

    // surfaceScalarField psi2_own =
    //     fvc::interpolate(psi2_, own_, "reconstruct(psi2)");

    // surfaceScalarField psi2_nei =
    //     fvc::interpolate(psi2_, nei_, "reconstruct(psi2)");

    // Wp_own_ = 1/(1.0 - (vF1face_*psi1_own + vF2face_*psi2_own)*ghf_);

    // Wp_nei_ = 1/(1.0 - (vF1face_*psi1_nei + vF2face_*psi2_nei)*ghf_);

    // rho0ghf_ = (vF1face_*rho01_ + vF2face_*rho02_)*ghf_;
}
