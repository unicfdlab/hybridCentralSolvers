{
    rho_pos     = fvc::interpolate(rho, pos, "reconstruct(rho)");
    rho_neg     = fvc::interpolate(rho, neg, "reconstruct(rho)");
    
    psi_pos     = fvc::interpolate(psi, pos, "reconstruct(psi)");
    psi_neg     = fvc::interpolate(psi, neg, "reconstruct(psi)");
    
    rhoHat_pos  = fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)");
    rhoHat_neg  = fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)");
    
    rhoU_pos    = fvc::interpolate(rho*HbyA, pos, "reconstruct(U)");
    rhoU_neg    = fvc::interpolate(rho*HbyA, neg, "reconstruct(U)");
    
    phiv_pos    = ((rhoU_pos / rho_pos) & mesh.Sf()) - meshPhi;
    phiv_neg    = ((rhoU_neg / rho_neg) & mesh.Sf()) - meshPhi;
    
    c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
    cSf_pos     = fvc::interpolate(c, pos, "reconstruct(c)")*mesh.magSf();
    cSf_neg     = fvc::interpolate(c, neg, "reconstruct(c)")*mesh.magSf();
    
    ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
    am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);
    
    a_pos       = ap/(ap - am);
    aSf         = am*a_pos;
    a_neg       = 1.0 - a_pos;
    
    aphiv_pos   = a_pos*phiv_pos - aSf;
    aphiv_neg   = a_neg*phiv_neg + aSf;
    
    phid_pos = aphiv_pos * psi_pos;
    phid_neg = aphiv_neg * psi_neg;
    
    Dp_pos = a_pos*fvc::interpolate(rho*rAU, pos, "reconstruct(Dp)");
    Dp_neg = a_neg*fvc::interpolate(rho*rAU, neg, "reconstruct(Dp)");
    
}

