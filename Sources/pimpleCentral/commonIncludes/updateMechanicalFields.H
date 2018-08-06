    K = 0.5*magSqr(U);
    
    //add convection part
    if (fluxScheme == "Tadmor")
    {
        surfaceVectorField U_own
        (
            "U_own",
            fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own
        );
        
        surfaceVectorField U_nei
        (
            "U_nei",
            fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei
        );

        KChange = 
            fvc::div(phi_own*0.5*magSqr(U_own))
            +
            fvc::div(phi_nei*0.5*magSqr(U_nei));
    }
    else if (fluxScheme == "Kurganov")
    {
        KChange = 
                fvc::div(phi_own,K)
                +
                fvc::div(phi_nei,K);
    }
    else
    {
    }
    KSource = 
            fvc::div
            (
                (
                    linearInterpolate
                    (
                        (-turbulence->devRhoReff()) & U
                    ) & mesh.Sf()
                )()
            );
    
    //add time derivative part
    if (LTS)
    {
        KChange +=
            (rho*K - rho.oldTime()*K.oldTime())*trDeltaT.ref();
        KSource +=
            (p - p.oldTime())*trDeltaT.ref();
    }
    else
    {
        KChange +=
            (rho*K - rho.oldTime()*K.oldTime())/runTime.deltaT();
        KSource +=
            (p - p.oldTime())/runTime.deltaT();
    }


