surfaceScalarField rhoHat_pos
(
    "rhoHat_pos",
    fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)")
);

surfaceScalarField rhoHat_neg
(
    "rhoHat_neg",
    fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)")
);

surfaceScalarField liqDiffusionFlux
(
    "liqDiffusionFlux",
    fvc::interpolate(turbulence->muEff()) * ( (fvc::snGrad(YLiq)) * mesh.magSf()) * 0.0
);

surfaceScalarField gasDiffusionFlux
(
    "gasDiffusionFlux",
    liqDiffusionFlux
);

surfaceScalarField phiRhoHat_pos
(
    "phiRhoHat_pos",
    phi_pos * 0.0
);
surfaceScalarField phiRhoHat_neg
(
    "phiRhoHat_neg",
    phiRhoHat_pos
);

surfaceScalarField YbarLiq_pos
(
    "YbarLiq_pos",
    fvc::interpolate(YbarLiq, pos, "reconstruct(YbarLiq)")
);

surfaceScalarField YbarLiq_neg
(
    "YbarLiq_neg",
    fvc::interpolate(YbarLiq, neg, "reconstruct(YbarLiq)")
);


