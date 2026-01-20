if (pimple.corrPISO() > 1)
{
    thermo.correct();
}

volScalarField rhoT = thermo.rho();
rhoHat = rhoT - psi*p;

rAU = 1.0 / UEqn.A();
HbyA = UEqn.H() * rAU;
HbyA.boundaryFieldRef() == U.boundaryField();

UpdateSonicSpeedOnFaces
(
    thermo,
    own,
    nei,
    p_own,
    p_nei,
    Rmf,
    cf_own,
    cf_nei
);

CalculateKnpWeights
(
    noFluxWeighting,
    fluxScheme,
    v_zero,
    thermo,
    mesh,
    rho, // Here we weight with rho
    U,
    c,
    psi,
    own,
    nei,
    mesh_phi,
    rho_own,
    rho_nei,
    phiv_own,
    phiv_nei,
    cf_own,
    cf_nei,
    cSf_own,
    cSf_nei,
    ap,
    am,
    alpha_own,
    alpha_nei,
    aSf
);

RhoUpdateCentralFields
(
    noFluxWeighting,
    MRF,
    mesh,
    rho,
    rhoHat,
    HbyA,
    rAU,
    psi,     // 1/(Rm*T)
    T,
    own,
    nei,
    mesh_phi,
    kappa,
    alpha_own,
    alpha_nei,
    rho_own,
    rho_nei,
    aSf,
    rhoU_own,
    rhoU_nei,
    rhoHatf,
    rhoHat_own,
    rhoHat_nei,
    phiRhoHat_own,
    phiRhoHat_nei,
    phiv_own,
    phiv_nei,
    aphiv_own,
    aphiv_nei,
    phid_own,
    phid_nei,
    psi_own,
    psi_nei,
    psif,
    Dp_own,
    Dp_nei
);

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

p.relax();

p_own = fvc::interpolate(p, own, "reconstruct(p)");
p_nei = fvc::interpolate(p, nei, "reconstruct(p)");

gradp = fvc::div((alpha_own*p_own + alpha_nei*p_nei)*mesh.Sf());
#include "limitPressureGradient.H"

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max(U): " << max(U).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

#include "massEqn.H"

rhoT = rhoHat + psi*p;

#include "centralContinuityErrs.H"

rho = rhoT;

rho_own = fvc::interpolate(rho, own, "reconstruct(rho)");
rho_nei = fvc::interpolate(rho, nei, "reconstruct(rho)");
correctCentralACMIInterpolation(rho_nei);
