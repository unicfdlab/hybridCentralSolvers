SuSp.set
(
    iCmpt,
    new fvScalarMatrix
    (
        reaction->R(Yi)
        +
        fvOptions(rho,Yi)
        +
        parcels.SYi(iCmpt, Yi)
    )
);
