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

lambdaCoeffs.primitiveFieldRef() = 1.0;

{
    word YiScheme ("div(phi," + YLiq.name() +")");
    
    surfaceScalarField rhoPhiYCorr
    (
        "rhoPhiYLiqCorr",
        fvc::flux
        (
            phi,
            YLiq,
            YiScheme
        )
    );
    
    surfaceScalarField mDiffC
    (
        "mDiffC",
        fvc::interpolate(turbulence->muEff()*0.0)
    );
    
    Info << "Creating matrix for fvOptions" << endl;
    
    #warning "change call to fvOptions to compressible"
    
    fvScalarMatrix SuSp
    (
        fvOptions(rho,YLiq)
    );
    
    mulesWithDiffusionImplicitLimiter
    (
        rho,
        YLiq,
        phi,
        allFacesLambda,
        rhoPhiYCorr,
        liqDiffusionFlux,
        mDiffC,
        SuSp
    );

    upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix YLiqEqn
    (
        fvm::ddt(rho,YLiq)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, YLiq)
        +
        fvc::div(lambdaCoeffs * rhoPhiYCorr)
        ==
        SuSp
    );
    
    YLiqEqn.relax();
    
    YLiqEqn.solve(mesh.solver(YLiq.name()));
    
    surfaceScalarField centralSurfaceFlux
    (
        "centralSurfaceFlux",
        fvc::flux
        (
            phi_own,
            YLiq,
            YiScheme
        )
        +
        fvc::flux
        (
            phi_nei,
            YLiq,
            YiScheme
        )
    );

    upwind<scalar> UDsPos(mesh, phi_own);
    upwind<scalar> UDsNeg(mesh, phi_nei);
    
    surfaceScalarField centralSurfaceFluxBD
    (
        fv::gaussConvectionScheme<scalar>(mesh, phi_own, UDsPos).flux(phi_own, YLiq)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi_nei, UDsNeg).flux(phi_nei, YLiq)
    );
    
    //account for diffusion in enthalpy equation
    centralSurfaceFlux =
        centralSurfaceFluxBD + lambdaCoeffs*(centralSurfaceFlux - centralSurfaceFluxBD);

    liqDiffusionFlux = -
        (
            centralSurfaceFlux - (YLiqEqn.flux() + lambdaCoeffs * rhoPhiYCorr)
        );
    
    gasDiffusionFlux = -liqDiffusionFlux;

    Info<< "Mass fraction min/max = " << min(YLiq).value() << "/"
        << max(YLiq).value() << endl;
    
    YLiq = max(YLiq, 0.0);
    YLiq = min(YLiq, 1.0);
    YGas = 1.0 - YLiq;
    
}
