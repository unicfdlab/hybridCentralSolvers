{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
    );
    
    rhoEqn.solve();
    
    Info<< "rho max/min : " 
        << gMax(rho.internalField())
        << " / " 
        << gMin(rho.internalField()) << endl;
}
