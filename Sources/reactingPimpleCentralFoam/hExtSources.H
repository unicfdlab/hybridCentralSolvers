fvScalarMatrix hExtSources
(
    reaction->Qdot()
    +
    fvOptions(rho,h)
);