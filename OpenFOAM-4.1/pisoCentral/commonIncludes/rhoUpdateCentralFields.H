{
    psi_own       = fvc::interpolate(psi, own, "reconstruct(psi)");
    psi_nei       = fvc::interpolate(psi, nei, "reconstruct(psi)");
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);
    
    rhoHat_own    = fvc::interpolate(rhoHat, own, "reconstruct(rhoHat)");
    rhoHat_nei    = fvc::interpolate(rhoHat, nei, "reconstruct(rhoHat)");
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(rhoHat_nei);
    
    rhoU_own      = fvc::interpolate(rho*HbyA, own, "reconstruct(U)");
    rhoU_nei      = fvc::interpolate(rho*HbyA, nei, "reconstruct(U)");
    
    phiv_own      = (rhoU_own / rho_own) & mesh.Sf();
    phiv_nei      = (rhoU_nei / rho_nei) & mesh.Sf();
    if (mesh.changing())
    {
        phiv_own   -= mesh.phi();
        phiv_nei   -= mesh.phi();
    }
    
    MRF.makeRelative(phiv_own);
    MRF.makeRelative(phiv_nei);
    
    aphiv_own     = alpha_own*phiv_own - aSf;
    aphiv_nei     = alpha_nei*phiv_nei + aSf;
    
    phid_own      = aphiv_own * psi_own;
    phid_nei      = aphiv_nei * psi_nei;
    
    phiRhoHat_own = aphiv_own * rhoHat_own;
    phiRhoHat_nei = aphiv_nei * rhoHat_nei;
    
    Dp_own = alpha_own*fvc::interpolate(rho*rAU, own, "reconstruct(Dp)");
    Dp_nei = alpha_nei*fvc::interpolate(rho*rAU, nei, "reconstruct(Dp)");
}

