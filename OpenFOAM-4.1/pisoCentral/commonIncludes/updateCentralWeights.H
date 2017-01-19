    {
        rho_own     = fvc::interpolate(rho, own, "reconstruct(rho)");
        rho_nei     = fvc::interpolate(rho, nei, "reconstruct(rho)");
        correctCentralACMIInterpolation(rho_nei);
        
        phiv_own    = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own) & mesh.Sf();
        phiv_nei    = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        
        c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
        cf_own      = fvc::interpolate(c, own, "reconstruct(psi)");
        cf_nei      = fvc::interpolate(c, nei, "reconstruct(psi)");
        correctCentralACMIInterpolation(cf_nei);
        cSf_own     = cf_own * mesh.magSf();
        cSf_nei     = cf_nei * mesh.magSf();
        
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
