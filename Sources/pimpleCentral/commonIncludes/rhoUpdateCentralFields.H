{
    psif = linearInterpolate(psi);
    psi_own       = fvc::interpolate(psi, own, "reconstruct(psi)")*kappa
                      + (1.0 - kappa)*psif;
    psi_nei       = fvc::interpolate(psi, nei, "reconstruct(psi)")*kappa
                      + (1.0 - kappa)*psif;
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);
    
    rhoHatf = linearInterpolate(rhoHat);
    rhoHat_own    = fvc::interpolate(rhoHat, own, "reconstruct(rhoHat)")*kappa
                      + (1.0 - kappa)*rhoHatf;
    rhoHat_nei    = fvc::interpolate(rhoHat, nei, "reconstruct(rhoHat)")*kappa
                      + (1.0 - kappa)*rhoHatf;
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(rhoHat_nei);
    
    if (noFluxWeighting)
    {
        phiv_own    = (fvc::interpolate(HbyA, own, "reconstruct(U)")) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(HbyA, nei, "reconstruct(U)")) & mesh.Sf();
    }
    else
    {
        rhoU_own      = fvc::interpolate(rho*HbyA, own, "reconstruct(U)");
        rhoU_nei      = fvc::interpolate(rho*HbyA, nei, "reconstruct(U)");
        phiv_own      = (rhoU_own / rho_own) & mesh.Sf();
        phiv_nei      = (rhoU_nei / rho_nei) & mesh.Sf();
    }
    
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

