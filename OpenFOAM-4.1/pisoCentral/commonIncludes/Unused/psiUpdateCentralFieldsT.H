{
    rho_own     = fvc::interpolate(rho, own, "reconstruct(rho)");
    rho_nei     = fvc::interpolate(rho, nei, "reconstruct(rho)");
    
    psi_own     = fvc::interpolate(psi, own, "reconstruct(psi)");
    psi_nei     = fvc::interpolate(psi, nei, "reconstruct(psi)");
    
    psiU_own    = fvc::interpolate(psi*HbyA, own, "reconstruct(U)");
    psiU_nei    = fvc::interpolate(psi*HbyA, nei, "reconstruct(U)");
    
    phiv_own    = (psiU_own / psi_own) & mesh.Sf();
    phiv_nei    = (psiU_nei / psi_nei) & mesh.Sf();
    
    MRF.makeRelative(phi_own);
    MRF.makeRelative(phi_nei);
    
    c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
    cSf_own     = fvc::interpolate(c, own, "reconstruct(psi)")*mesh.magSf();
    cSf_nei     = fvc::interpolate(c, nei, "reconstruct(psi)")*mesh.magSf();
    
    ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
    am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);
    
    alpha_own   = 0.5;
    aSf         = -0.5*max(mag(ap), mag(am));
    alpha_nei   = 1.0 - alpha_own;
    
    aphiv_own   = alpha_own*phiv_own - aSf;
    aphiv_nei   = alpha_nei*phiv_nei + aSf;

    phid_own    = aphiv_own * psi_own;
    phid_nei    = aphiv_nei * psi_nei;
    
    Dp_own = alpha_own*fvc::interpolate(rho*rAU, own, "reconstruct(Dp)");
    Dp_nei = alpha_nei*fvc::interpolate(rho*rAU, nei, "reconstruct(Dp)");
}

