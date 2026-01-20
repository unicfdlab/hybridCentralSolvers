fvScalarMatrix hExtSources
(
    parcels.Sh(h)
    +
    reaction->Qdot()
    +
    fvOptions(rho,h)
);

