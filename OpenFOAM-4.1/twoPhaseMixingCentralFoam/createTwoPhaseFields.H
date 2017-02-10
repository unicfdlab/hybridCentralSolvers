surfaceScalarField mDCf
(
    "mDCf",
    phi*dimensionedScalar("zerom", dimless/dimLength, 0.0)
);

surfaceScalarField liqDiffusionFlux
(
    "liqDiffusionFlux",
    mDCf * fvc::snGrad(YLiq) * mesh.magSf() * 0.0
);

surfaceScalarField gasDiffusionFlux
(
    "gasDiffusionFlux",
    mDCf * fvc::snGrad(YGas) * mesh.magSf() * 0.0
);

