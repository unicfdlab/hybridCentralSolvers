    {
        fvScalarMatrix rhoEqn
        (
            fvm::ddt(rho) + fvc::div(phi)
            ==
            fvOptions(rho)
        );
        
        fvOptions.constrain(rhoEqn);
        
        rhoEqn.solve();
        
        fvOptions.correct(rho);
        
        Info<< "rho max/min : " << max(rho).value()
            << " / " << min(rho).value() << endl;
    }
    