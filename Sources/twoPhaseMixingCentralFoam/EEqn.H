{
    volScalarField& he = thermo.he();
    volScalarField alphahEff (turbulence->alphaEff());
    volScalarField Cp = thermo.Cp();

    upwind<scalar> UDsPos(mesh, phi_own);
    upwind<scalar> UDsNeg(mesh, phi_nei);
    upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix hConvection
    (
        fvm::ddt(rho,he)
        //+
        //fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, he)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi_own, UDsPos).fvmDiv(phi_own, he)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi_nei, UDsNeg).fvmDiv(phi_nei, he)
    );

    /*
     *
     * Explicit flux correction, bounded with MULES
     *
     */
    
    surfaceScalarField hFluxCorr
    (
        lambdaCoeffs
        *
        (
	    fvc::flux
	    (
		phi_own,
		he,
		"div(phi_own," + he.name() + ")"
	    )
	    +
	    fvc::flux
	    (
		phi_nei,
		he,
		"div(phi_nei," + he.name() + ")"
	    )
	    -
	    hConvection.flux()
        )
    );
    
    surfaceScalarField hCmptFlux
    (
        "hCmptFlux",
        fvc::flux
        (
            liqDiffusionFlux,
            thermo.thermoLiq().he(),
            "div(rhoi*Uri,hi)"
        )
        +
        fvc::flux
        (
            gasDiffusionFlux,
            thermo.thermoGas().he(),
            "div(rhoi*Uri,hi)"
        )
        +
        fvc::flux
        (
            fvc::snGrad(YLiq) * mesh.magSf(),
            thermo.thermoLiq().he() * alphahEff,
            "div(alphaEff*hi,grad(Yi))"
        )
        +
        fvc::flux
        (
            fvc::snGrad(YGas) * mesh.magSf(),
            thermo.thermoGas().he() * alphahEff,
            "div(alphaEff*hi,grad(Yi))"
        )
    );
    
    fvScalarMatrix EEqn
    (
        hConvection + fvc::div(hFluxCorr) + fvc::div(hCmptFlux) 
        +
        KChange//EkChange
        ==
        fvm::laplacian(alphahEff, he)
        +
        KSource
//        +
//        rho*(g & U)
//        +
//        fvOptions(rho,he)
    );
    
    fvOptions.constrain(EEqn);
    
    EEqn.solve();
    
    fvOptions.correct(he);
    
    thermo.correct();
    
    //surfaceScalarField eFlux = EEqn.flux() + hCmptFlux - hConvection.flux() - hFluxCorr;
    //
    //Info << "Max/min of eFlux: " << max(eFlux).value() << "/" << min(eFlux).value() << endl;
    //he.write();
    //alphahEff.rename("alphahEff");
    //alphahEff.write();
    
    //T.correctBoundaryConditions();
    
}
