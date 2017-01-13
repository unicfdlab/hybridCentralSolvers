fvScalarMatrix hEqn
(
    fvm::ddt(rho,h)
    +
    fvm::div(phi_own,h)
    +
    fvm::div(phi_nei,h)
    -
    fvm::laplacian(turbulence->alphaEff(), h)
    +
    KChange
    ==
    KSource
    +
    fvOptions(rho,h)
);

fvOptions.constrain(hEqn);

hEqn.solve();

fvOptions.correct(h);

thermo.correct();

Info << "max/min T: " << max(T).value() << "/" << min(T).value() << endl;
