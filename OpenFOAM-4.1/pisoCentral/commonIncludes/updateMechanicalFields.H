    K = 0.5*magSqr(U);
    
    {
        if (fluxScheme == "Tadmor")
        {
            surfaceVectorField U_own
            (
                "U_own",
                fvc::interpolate(rho*U, own, "reconstruct(rho*U)") / rho_own
            );
            
            surfaceVectorField U_nei
            (
                "U_nei",
                fvc::interpolate(rho*U, nei, "reconstruct(rho*U)") / rho_nei
            );

            KChange = 
                fvc::ddt(rho,K)
                +
                fvc::div(phi_own*0.5*magSqr(U_own))
                +
                fvc::div(phi_nei*0.5*magSqr(U_nei));
        }
        else //Kurganov
        {
            KChange = 
                fvc::ddt(rho,K)
                +
                fvc::div(phi_own,K)
                +
                fvc::div(phi_nei,K);
        }
    }
    
    KSource = fvc::ddt(p) + fvc::div
            (
                (
                    linearInterpolate
                    (
                        (-turbulence->devRhoReff()) & U
                    ) & mesh.Sf()
                )()
            );

    if (mesh.changing())
    {
        KSource -= fvc::div(fvc::meshPhi(rho,U),p);
    }

