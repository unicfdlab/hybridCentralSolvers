
volScalarField psiCpmCvByT = psi* (thermo.Cp() - thermo.Cv()) / T;

surfaceScalarField pressureGradEnergyFluxContribution
(
    "pressureGradEnergyFluxContribution",
    linearInterpolate
    (
        (turbulence->alphaEff()/rho)*
        (1.0 - T*sqrt(psiCpmCvByT))*
        gradp
    ) & mesh.Sf()
);

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
    fvc::div(pressureGradEnergyFluxContribution)
    +
    KSource
    +
    fvOptions(rho,h)
);

hEqn.relax();

fvOptions.constrain(hEqn);

hEqn.solve();

fvOptions.correct(h);

thermo.correct();

Info << "max/min T: " << max(T).value() << "/" << min(T).value() << endl;
