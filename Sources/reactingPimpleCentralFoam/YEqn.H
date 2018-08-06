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
    reaction->correct();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);
    const label nComponents = Y.size();

    PtrList<fvScalarMatrix> SuSp(nComponents);
    
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
            mDCf[iCmpt] = fvc::interpolate(turbulence->muEff());
        }
        else
        {
            inertIndex = iCmpt;
            forAll(diffusiveFlux[iCmpt].internalField(), iFace)
            {
                diffusiveFlux[iCmpt].primitiveFieldRef()[iFace] = 0.0;
            }
            forAll(diffusiveFlux[iCmpt].boundaryFieldRef(), iPatch)
            {
                forAll(diffusiveFlux[iCmpt].boundaryFieldRef()[iPatch], iFace)
                {
                    diffusiveFlux[iCmpt].boundaryFieldRef()[iPatch][iFace] = 0.0;
                }
            }
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
                fvm::ddt(rho,Yi)
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
            
            YiEqn.relax();
            
            fvOptions.constrain(YiEqn);
            
            YiEqn.solve(mesh.solver("Yi"));
            
            fvOptions.correct(Yi);
            
            Yi.max(0.0);
            Yt += Yi;
        
            phiY[iCmpt] = YConvection.flux() + lambdaCoeffs * rhoPhiYCorr;
            
            diffusiveFlux[iCmpt] = YDiffusion.flux();
            diffusiveFlux[inertIndex] -= diffusiveFlux[iCmpt];
            
            Info << "max/min " << Yi.name() << ": " << max(Yi).value() << "/" << min(Yi).value() << endl;
        }
    }
    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].correctBoundaryConditions();
    Y[inertIndex].max(0.0);
}
