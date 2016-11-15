    {
        phiv_own    = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        
        //phiv_own    = (fvc::interpolate(rho*(HbyA - rAU*gradp), own, "reconstruct(U)") / rho_own) & mesh.Sf();
        //phiv_nei    = (fvc::interpolate(rho*(HbyA - rAU*gradp), nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        
        c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
        cSf_own     = fvc::interpolate(c, own, "reconstruct(psi)")*mesh.magSf();
        cSf_nei     = fvc::interpolate(c, nei, "reconstruct(psi)")*mesh.magSf();
        
        ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
        am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);
        
        if (fluxScheme == "Kurganov")
        {
            alpha_own   = ap/(ap - am);
            aSf         = am*alpha_own;
            alpha_nei   = 1.0 - alpha_own;
        }
        else if (fluxScheme == "Tadmor")
        {
            aSf         = -0.5*max(mag(ap), mag(am));
        }
        else
        {
        
        }
    }
