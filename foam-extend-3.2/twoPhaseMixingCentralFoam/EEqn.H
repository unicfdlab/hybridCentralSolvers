{
    volScalarField& he = thermo.h();
    volScalarField alphahEff (turbulence->alphaEff());
    volScalarField Cp = thermo.Cp();

    upwind<scalar> UDsPos(mesh, phi_pos);
    upwind<scalar> UDsNeg(mesh, phi_neg);
    //upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix hConvection
    (
        fvm::ddt(rho,he)
        //+
        //fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, he)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi_pos, UDsPos).fvmDiv(phi_pos, he)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi_neg, UDsNeg).fvmDiv(phi_neg, he)
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
		phi_pos,
		he,
		"div(phiPos," + he.name() + ")"
	    )
	    +
	    fvc::flux
	    (
		phi_neg,
		he,
		"div(phiNeg," + he.name() + ")"
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
        KChange
        ==
        fvm::laplacian(alphahEff, he)
        +
        dpdt
//        +
//        rho*(g & U)
    );
    
    EEqn.solve();
    
    thermo.correct();
}
