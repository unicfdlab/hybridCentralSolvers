    {
        if (noFluxWeighting)
        {
            phiv_own    = (fvc::interpolate(U, own, "reconstruct(U)")) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(U, nei, "reconstruct(U)")) & mesh.Sf();
        }
        else
        {
            rho_own     = fvc::interpolate(rho, own, "reconstruct(rho)");
            rho_nei     = fvc::interpolate(rho, nei, "reconstruct(rho)");
            correctCentralACMIInterpolation(rho_nei);
            phiv_own    = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own) & mesh.Sf();
            phiv_nei    = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei) & mesh.Sf();
        }
        
        if (mesh.changing())
        {
            phiv_own   -= mesh.phi();
            phiv_nei   -= mesh.phi();
        }
        
        c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
        cf_own      = fvc::interpolate(c, own, "reconstruct(psi)");
        cf_nei      = fvc::interpolate(c, nei, "reconstruct(psi)");
        correctCentralACMIInterpolation(cf_nei);
        cSf_own     = cf_own * mesh.magSf();
        cSf_own.setOriented(true);
        cSf_nei     = cf_nei * mesh.magSf();
        cSf_nei.setOriented(true);
        
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
            FatalErrorIn
            (
                    "updateCentralWeights.H"
            )   << "Unknown central flux scheme " << fluxScheme << nl << nl
                << "Valid flux schemes are:" << nl
                << "1) Kurganov" << nl
                << "2) Tamodr" << endl;
        }
    }
