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

surfaceScalarField hLambdaCoeffs ("hLambdaCoeffs", own);
surfaceScalarField maxDeltaY ("maxDeltaY", own*0.0);

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
            mDCf[iCmpt] = linearInterpolate(turbulence->muEff());
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
                YConvection + YDiffusion// +  fvc::div(lambdaCoeffs * rhoPhiYCorr)
            );

            //YiEqn.relax();

            //fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            //fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;

            //phiY[iCmpt] = YConvection.flux() + lambdaCoeffs * rhoPhiYCorr;

            diffusiveFlux[iCmpt] = YDiffusion.flux();
            diffusiveFlux[inertIndex] -= diffusiveFlux[iCmpt];

            Info << "max/min " << Yi.name() << ": " << max(Yi).value() << "/" << min(Yi).value() << endl;
        }
    }
    Y[inertIndex].primitiveFieldRef() = (scalar(1.0) - Yt.primitiveField());
    Y[inertIndex].correctBoundaryConditions();
    Y[inertIndex].max(0.0);

    label owner = -1;
    label neighbour = -1;
    scalar deltaY = 0.0;
    forAll(Y, iCmpt)
    {
        forAll(maxDeltaY.primitiveField(), iFace)
        {
            if (mesh.isInternalFace(iFace))
            {
                owner = mesh.owner()[iFace];
                neighbour = mesh.neighbour()[iFace];
                deltaY = mag(Y[iCmpt][owner] - Y[iCmpt][neighbour]);
                if (deltaY > maxDeltaY.primitiveField()[iFace])
                {
                     maxDeltaY.primitiveFieldRef()[iFace] = deltaY;
                }
            }
        }
        forAll(maxDeltaY.boundaryField(), iPatch)
        {
            if (maxDeltaY.boundaryField()[iPatch].coupled())
            {
                scalarField intF  (maxDeltaY.boundaryField()[iPatch]);
                scalarField intY  (Y[iCmpt].boundaryField()[iPatch]);
                scalarField neiY  (Y[iCmpt].boundaryField()[iPatch].patchNeighbourField());
                forAll(intF, iFace)
                {
                    deltaY = mag(intY[iFace] - neiY[iFace]);
                    if (deltaY > intF[iFace])
                    {
                        intF[iFace] = deltaY;
                    }
                }
                maxDeltaY.boundaryFieldRef()[iPatch].operator = (intF);
            }
        }
    }
    forAll(maxDeltaY.primitiveField(), iFace)
    {
        if (maxDeltaY.primitiveField()[iFace] >= 0.05)
        {
            hLambdaCoeffs.primitiveFieldRef()[iFace] = 
                lambdaCoeffs.primitiveField()[iFace];
        }
    }
    forAll(maxDeltaY.boundaryField(), iPatch)
    {
        if (maxDeltaY.boundaryField()[iPatch].coupled())
        {
            scalarField intF = maxDeltaY.boundaryField()[iPatch];
            scalarField intH = hLambdaCoeffs.boundaryField()[iPatch];
            const scalarField& intL = lambdaCoeffs.boundaryField()[iPatch];
            forAll(intF, iFace)
            {
                if (intF[iFace] > 0.05)
                {
                    intH[iFace] = intL[iFace];
                }
            }
            hLambdaCoeffs.boundaryFieldRef()[iPatch].operator = (intH);
        }
    }
}
