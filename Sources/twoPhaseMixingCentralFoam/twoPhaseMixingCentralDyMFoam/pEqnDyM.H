
if (pimple.corrPISO() > 1 )
{
    thermo.correct();
}

volScalarField rhoT = thermo.rho();
rhoHat = rhoT - psi*p;
volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());

HbyA = UEqn.H() * rAU;

/*
 *
 * Hybrid FVM p-U coupling procedure
 *
 */
{
    /*
    rho_own       = fvc::interpolate(rho, own, "reconstruct(rho)");
    rho_nei       = fvc::interpolate(rho, nei, "reconstruct(rho)");
    correctCentralACMIInterpolation(rho_nei);
    
    psi_own       = fvc::interpolate(psi, own, "reconstruct(psi)");
    psi_nei       = fvc::interpolate(psi, nei, "reconstruct(psi)");
    correctCentralACMIInterpolation(psi_nei);
    
    rhoHat_own    = fvc::interpolate(rhoHat, own, "reconstruct(rhoHat)");
    rhoHat_nei    = fvc::interpolate(rhoHat, nei, "reconstruct(rhoHat)");
    correctCentralACMIInterpolation(rhoHat_nei);
    
    rhoU_own      = fvc::interpolate(rho*HbyA, own, "reconstruct(U)");
    rhoU_nei      = fvc::interpolate(rho*HbyA, nei, "reconstruct(U)");

    phiv_own      = (rhoU_own / rho_own) & mesh.Sf();
    phiv_nei      = (rhoU_nei / rho_nei) & mesh.Sf();
    
    c             = sqrt(thermo.Cp()/thermo.Cv() / psi);
    cf_own        = fvc::interpolate(c, own, "reconstruct(c)");
    cf_nei        = fvc::interpolate(c, nei, "reconstruct(c)");
    cSf_own       = cf_own*mesh.magSf();
    cSf_nei       = cf_nei*mesh.magSf();
    
    phiv_own -= meshPhi;
    phiv_nei -= meshPhi;
    
    MRF.makeRelative(phiv_own);
    MRF.makeRelative(phiv_nei);
    
    ap            = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
    am            = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);
    
    alpha_own     = ap/(ap - am);
    aSf           = am*alpha_own;
    alpha_nei     = 1.0 - alpha_own;
    
    aphiv_own     = alpha_own * phiv_own - aSf;
    aphiv_nei     = alpha_nei * phiv_nei + aSf;
    
    phid_own      = aphiv_own * psi_own;
    phid_nei      = aphiv_nei * psi_nei;
    
    phiRhoHat_own = aphiv_own * rhoHat_own;
    phiRhoHat_nei = aphiv_nei * rhoHat_nei;
    
    Dp_own = alpha_own*fvc::interpolate(rho*rAU, own, "reconstruct(Dp)");

    Dp_nei = alpha_nei*fvc::interpolate(rho*rAU, nei, "reconstruct(Dp)");
    */
    
    //#include "updateCentralWeightsDyM.H"
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
            phi_own = pEqn_own.flux() + phiRhoHat_own;//Relative
            phi_nei = pEqn_nei.flux() + phiRhoHat_nei;//Relative
            phi = phi_own + phi_nei;                  //Relative
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

p = max(p,pMin);
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

