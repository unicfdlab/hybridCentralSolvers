    /*
     *
     * Multicomponent fields
     *
     */
    // Convective and Diffusive mass fluxes
    PtrList<surfaceScalarField> mDCf(Y.size());
    PtrList<surfaceScalarField> diffusiveFlux(Y.size());
    PtrList<surfaceScalarField> phiY(Y.size());
    for (label cmptIdx=0; cmptIdx < Y.size(); cmptIdx++)
    {
        mDCf.set
        (
            cmptIdx,
            new surfaceScalarField
            (
                "mDCf",
                fvc::interpolate(turbulence->muEff())
            )
        );
        
        diffusiveFlux.set
        (
            cmptIdx,
            new surfaceScalarField
            (
                "phiD",
                mDCf[cmptIdx] * fvc::snGrad(Y[cmptIdx]) * mesh.magSf() * 0.0
            )
        );
        
        if (Y[cmptIdx].name() == inertSpecie)
        {
            phiY.set
            (
                cmptIdx,
                new surfaceScalarField
                (
                    "rhoiPhi",
                    phi * 0.0
                )
            );
        }
    }
//
// END-OF-FILE
//