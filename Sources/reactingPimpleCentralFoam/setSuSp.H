SuSp.set
(
    iCmpt,
    new fvScalarMatrix
    (
        reaction->R(Yi)
        +
        fvOptions(rho,Yi)
    )
);
