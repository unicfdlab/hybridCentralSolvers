if (corr > 1)
{
    thermo.correct();
}

volScalarField rhoT = thermo.rho();
rhoHat = rhoT - psi*p;

volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());

HbyA = UEqn.H() * rAU;
{
    #include "rhoUpdateCentralFields.H"

    YbarLiq_pos = fvc::interpolate(YbarLiq, pos, "reconstruct(YbarLiq)");
    YbarLiq_neg = fvc::interpolate(YbarLiq, neg, "reconstruct(YbarLiq)");

    phiRhoHat_pos = aphiv_pos * rhoHat_pos;
    phiRhoHat_neg = aphiv_neg * rhoHat_neg;

    //update anti-diffusion fluxes in the interface regions
    {
        phiRhoHat_pos -=
            aSf * (YbarLiq_neg - YbarLiq_pos) * rhoHat_pos * Cid;
        
        phiRhoHat_neg -=
            aSf * (YbarLiq_neg - YbarLiq_pos) * rhoHat_neg * Cid;
        
    }

    for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn_pos
        (
            fvm::div(phid_pos,p) + fvc::div(phiRhoHat_pos) - fvm::laplacian(Dp_pos,p)
        );
    
        fvScalarMatrix pEqn_neg
        (
            fvm::div(phid_neg,p) + fvc::div(phiRhoHat_neg) - fvm::laplacian(Dp_neg,p)
        );
    
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
            +
            fvc::ddt(rhoHat)
            +
            pEqn_pos
            +
            pEqn_neg
        );
    
        if
        (
            oCorr == nOuterCorr-1
            && corr == nCorr-1
            && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve(mesh.solutionDict().solver("pFinal"));
        }
        else
        {
            pEqn.solve();
        }
    
        if (nonOrth == nNonOrthCorr)
        {
            phi_pos = pEqn_pos.flux() + phiRhoHat_pos;
            phi_neg = pEqn_neg.flux() + phiRhoHat_neg;
            phi = phi_pos + phi_neg;
        }
    }

    p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
    p_neg = fvc::interpolate(p, neg, "reconstruct(p)");
    rho_pos = p_pos*psi_pos + rhoHat_pos;
    rho_neg = p_neg*psi_neg + rhoHat_neg;

    gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());

    #include "limitPressureGradient.H"
}

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();

K = 0.5*magSqr(U);

Info << "max(mag(U)): " << max(mag(U)).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

#include "MixtureRhoEqn.H"

p = max(p, pMin);
p.correctBoundaryConditions();
rhoT = rhoHat + psi*p; //1.0 / (YLiq / rhoLiq + YGas / rhoGas);

#include "centralContinuityErrs.H"

rho = rhoT;


