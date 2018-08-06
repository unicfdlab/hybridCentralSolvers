volScalarField rhoHat
(
    IOobject
    (
        "rhoHat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho  - psi*p
);

surfaceScalarField rhoHatf
(
    "rhoHatf",
    linearInterpolate(rhoHat)
);

surfaceScalarField rhoHat_own
(
    "rhoHat_own",
    fvc::interpolate(rhoHat, own, "reconstruct(rhoHat)")
);

surfaceScalarField rhoHat_nei
(
    "rhoHat_nei",
    fvc::interpolate(rhoHat, nei, "reconstruct(rhoHat)")
);

surfaceScalarField phiRhoHat_own
(
    "phiRhoHat_own",
    rhoHat_own * aphiv_own
);

surfaceScalarField phiRhoHat_nei
(
    "phiRhoHat_nei",
    rhoHat_nei * aphiv_nei
);

surfaceVectorField rhoU_own
(
    "rhoU_own",
    fvc::interpolate(rho*U, own, "reconstruct(U)")
);

surfaceVectorField rhoU_nei
(
    "rhoU_nei",
    fvc::interpolate(rho*U, nei, "reconstruct(U)")
);
