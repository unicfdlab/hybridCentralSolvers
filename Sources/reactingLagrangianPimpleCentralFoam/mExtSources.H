fvScalarMatrix mExtSources
(
    parcels.Srho()
    +
    fvOptions(psi, p, rho.name())
);
