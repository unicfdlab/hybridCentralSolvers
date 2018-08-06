Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(fluidMesh)
);
psiThermo& thermo = pThermo();
thermo.correct();

volScalarField& p = thermo.p();
volScalarField& h = thermo.he();
volScalarField& psi = const_cast<volScalarField&>(thermo.psi());
const volScalarField& T = thermo.T();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        fluidMesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    fluidMesh
);

volVectorField HbyA ("HbyA", U);


Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K
(
    "K",
    0.5*magSqr(U)
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        fluidMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

dimensionedScalar v_one("v_one", pow(dimTime,-1.0), 1.0);
volScalarField rAU ("rAU", 1.0/rho/v_one);
rAU.primitiveFieldRef() *= runTime.deltaTValue();
rAU.primitiveFieldRef() /= mesh.V().field();


volScalarField KChange
(
    "KChange",
    rho*K / runTime.deltaT() * 0.0
);

volScalarField KSource
(
    "KSource",
    p / runTime.deltaT() * 0.0
);

volScalarField c
(
    "cSound",
    sqrt(thermo.Cp()/thermo.Cv() / psi)
);

volVectorField gradp
(
    "gradp",
    fvc::grad(p)
);
volVectorField limGradP 
(
    "limGradP",
    fvc::grad(p)
);

DynamicList<label> badQualityCells(0);

Switch updateEnergyInPISO (false);
word fluxScheme (word::null);
Switch noFluxWeighting(false);

#include "readAdditionalPimpleControl.H"
#include "createCommonCentralFields.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

#include "createFvOptions.H"
#include "createMRF.H"
#include "initContinuityErrs.H"
#include "readCourantType.H"

#include "markBadQualityCells.H"

#include "psiUpdateCentralFields.H"
#include "updateKappa.H"

#include "updateCentralWeights.H"
#include "createCentralCourantNo.H"
