// Create limiter field for mass and energy fluxes
scalarField allFacesLambda(mesh.nFaces(), 1.0);
slicedSurfaceScalarField lambdaCoeffs
(
    IOobject
    (
        "lambdaCoeffs",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    ),
    mesh,
    dimless,
    allFacesLambda,
    false   // Use slices for the couples
);

upwind<scalar> UDsOwn(mesh, phi_own);
upwind<scalar> UDsNei(mesh, phi_nei);

{
    /*
     *
     * Create fluxes and sources
     *
     */

    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    Yt.rename("sumNonInertY");

    PtrList<fvScalarMatrix> SuSp(Y.size());
    forAll(Y, iCmpt)
    {
        volScalarField& Yi = Y[iCmpt];
        if (Yi.name() != inertSpecie)
        {
            word YiOwnScheme ("div(phi_own,Yi)");
            word YiNeiScheme ("div(phi_nei,Yi)");

            phiY.set
            (
                iCmpt,
                fvc::flux
                (
                    phi_own,
                    Yi,
                    YiOwnScheme
                )
                +
                fvc::flux
                (
                    phi_nei,
                    Yi,
                    YiNeiScheme
                )
            );
            #include "setSuSp.H"
            mDCf[iCmpt] = linearInterpolate(turbulence->muEff() / ScNumbers[iCmpt]);
        }
        else
        {
            inertIndex = iCmpt;
        }
    }


    /*
     *
     * Calculate the limiter
     *
     */
    surfaceScalarField phi_nei_zero = phi_nei*0.0;
    forAll(Y, iCmpt)
    {
        volScalarField& Yi = Y[iCmpt];
        if ( Yi.name() != inertSpecie )
        {
            surfaceScalarField& rhoPhiYCorr = phiY[iCmpt];
            mulesWithDiffusionImplicitLimiter
            (
                rho,
                Yi,
                phi_own,
                phi_nei,
                allFacesLambda,
                rhoPhiYCorr,
                diffusiveFlux[iCmpt],
                mDCf[iCmpt],
                SuSp[iCmpt]
            );
        }
    }

    /*
     *
     * Solve for components
     *
     */
    volScalarField MuInv = (Rm / Foam::constant::physicoChemical::R)*0.0;
    diffusiveFlux[inertIndex].boundaryFieldRef() = 0.0;
    diffusiveFlux[inertIndex].primitiveFieldRef() = 0.0;
    forAll(Y, iCmpt)
    {
        volScalarField& Yi = Y[iCmpt];
        if (Yi.name() != inertSpecie)
        {
            surfaceScalarField& rhoPhiYCorr = phiY[iCmpt];
            fvScalarMatrix YConvection
            (
                fvm::ddt(rho, Yi)
                +
                fv::gaussConvectionScheme<scalar>(mesh, phi_own, UDsOwn).fvmDiv(phi_own, Yi)
                +
                fv::gaussConvectionScheme<scalar>(mesh, phi_nei, UDsNei).fvmDiv(phi_nei, Yi)
                ==
                SuSp[iCmpt]
            );

            fvScalarMatrix YDiffusion
            (
                - fvm::laplacian(mDCf[iCmpt], Yi)
            );

            fvScalarMatrix YiEqn
            (
                YConvection + YDiffusion +  fvc::div(lambdaCoeffs * rhoPhiYCorr)
            );

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;

            diffusiveFlux[iCmpt] = YDiffusion.flux();
            diffusiveFlux[inertIndex] -= diffusiveFlux[iCmpt];

            dimensionedScalar MuCmptInv ("MuCmpt", MuInv.dimensions(), 1000.0 /thermo.composition().W(iCmpt));
            MuInv += (Yi * MuCmptInv);
        }
    }
    Y[inertIndex].primitiveFieldRef() = (scalar(1.0) - Yt.primitiveField());
    Y[inertIndex].correctBoundaryConditions();
    Y[inertIndex].max(0.0);
    dimensionedScalar MuCmptInv ("MuCmpt", MuInv.dimensions(), 1000.0 / thermo.composition().W(inertIndex));
    MuInv += (Y[inertIndex] * MuCmptInv);
    Rm = Foam::constant::physicoChemical::R * MuInv;
}
