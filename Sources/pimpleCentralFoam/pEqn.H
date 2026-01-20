rAU = 1.0 / UEqn.A();
HbyA = UEqn.H() * rAU;
HbyA.boundaryFieldRef() == U.boundaryField();

W = rho*Rm;

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
    W, // use p/T instead of rho
    U,
    c,
    psi,
    own,
    nei,
    mesh_phi,
    W_own,
    W_nei,
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

PsiUpdateCentralFields
(
    noFluxWeighting,
    MRF,
    mesh,
    W, // p/T
    HbyA,
    rAU,
    psi,     // 1/(Rm*T)
    T,
    Rm,
    own,
    nei,
    mesh_phi,
    kappa,
    alpha_own,
    alpha_nei,
    aSf,
    Rmf,
    W_own,
    W_nei,
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
        fvm::div(phid_own,p) - fvm::laplacian(Dp_own,p)
    );

    fvScalarMatrix pEqn_nei
    (
        fvm::div(phid_nei,p) - fvm::laplacian(Dp_nei,p)
    );

    fvScalarMatrix pEqn
    (
        fvm::ddt(psi,p)
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
        phi_own = pEqn_own.flux();
        phi_nei = pEqn_nei.flux();
        phi = phi_own + phi_nei;
    }
}

p.relax();
p_own = fvc::interpolate(p, own, "reconstruct(p)");
p_nei = fvc::interpolate(p, nei, "reconstruct(p)");

rho_own = p_own*psi_own;
rho_nei = p_nei*psi_nei;
correctCentralACMIInterpolation(rho_nei);

gradp = fvc::div((alpha_own*p_own + alpha_nei*p_nei)*mesh.Sf());
#include "limitPressureGradient.H"

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max(U): " << max(U).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

#include "massEqn.H"

volScalarField rhoT (thermo.rho());

#include "centralContinuityErrs.H"

rho = rhoT;
