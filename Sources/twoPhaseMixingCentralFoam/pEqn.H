if (pimple.corrPISO() > 1)
{
    thermo.correct();
}

volScalarField rhoT = thermo.rho();
rhoHat = rhoT - psi*p;
volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());

HbyA = UEqn.H() * rAU;

{
    #include "updateCentralWeights.H"
    #include "rhoUpdateCentralFields.H"

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn_own
        (
            fvm::div(phid_own,p) + fvc::div(phiRhoHat_own) - fvm::laplacian(Dp_own,p)
        );
    
        fvScalarMatrix pEqn_nei
        (
            fvm::div(phid_nei,p) + fvc::div(phiRhoHat_nei) - fvm::laplacian(Dp_nei,p)
        );
    
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
            +
            fvc::ddt(rhoHat)
            +
            pEqn_own
            +
            pEqn_nei
            ==
            fvOptions(psi,p,rho.name())
        );
    
        fvOptions.constrain(pEqn);
    
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
        if (pimple.finalNonOrthogonalIter())
        {
            phi_own = pEqn_own.flux() + phiRhoHat_own;
            phi_nei = pEqn_nei.flux() + phiRhoHat_nei;
            phi = phi_own + phi_nei;                  
        }
    }
    
    p_own = fvc::interpolate(p, own, "reconstruct(p)");
    p_nei = fvc::interpolate(p, nei, "reconstruct(p)");
    
    gradp = fvc::div((alpha_own*p_own + alpha_nei*p_nei)*mesh.Sf());

    #include "limitPressureGradient.H"
}

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max(mag(U)): " << max(mag(U)).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

p = max(p, pMin);
thermo.correctRealDensities();
rhoT = 1.0 / (YLiq / thermo.thermoLiq().rho() + YGas / thermo.thermoGas().rho());
#include "massEqn.H"

#include "centralContinuityErrs.H"
rho = rhoT;

rho_own = fvc::interpolate(rho, own, "reconstruct(rho)");
rho_nei = fvc::interpolate(rho, nei, "reconstruct(rho)");
correctCentralACMIInterpolation(rho_nei);

//
//END-OF-FILE
//

