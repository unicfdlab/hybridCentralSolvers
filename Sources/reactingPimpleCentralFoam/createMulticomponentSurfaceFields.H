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
            "mDCf" + Y[cmptIdx].name(),
            fvc::interpolate(turbulence->muEff())
        )
    );

    diffusiveFlux.set
    (
        cmptIdx,
        new surfaceScalarField
        (
            "phiD" + Y[cmptIdx].name(),
            mDCf[cmptIdx] * fvc::snGrad(Y[cmptIdx]) * mesh.magSf() * 0.0
        )
    );
}
//
// END-OF-FILE
//

