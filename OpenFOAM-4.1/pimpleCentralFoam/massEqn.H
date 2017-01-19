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
        
        Info<< "rho max/min : " << gMax(rho.primitiveField())
            << " / " << gMin(rho.primitiveField()) << endl;
    }
    