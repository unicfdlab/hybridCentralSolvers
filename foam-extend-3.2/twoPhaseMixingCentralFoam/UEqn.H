    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi_pos, U)
      + fvm::div(phi_neg, U)
      //+ fvm::div(liqDiffusionFlux, U)
      //+ fvm::div(gasDiffusionFlux, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
    );
    
    if (momentumPredictor)
    {
        solve
        (
            UEqn == -gradp
        );
        
        K = 0.5*magSqr(U);
        
        Info << "max(mag(U)): " << max(mag(U)).value() << endl;
    }
    
