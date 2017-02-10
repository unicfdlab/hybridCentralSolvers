    K = 0.5*magSqr(U);
    
    {
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
                (rho*K - rho.oldTime()*K.oldTime())/runTime.deltaT()
                +
                fvc::div(phi_own*0.5*magSqr(U_own))
                +
                fvc::div(phi_nei*0.5*magSqr(U_nei));
        }
        else //Kurganov
        {
            KChange = 
                (rho*K - rho.oldTime()*K.oldTime())/runTime.deltaT()
                +
                fvc::div(phi_own,K)
                +
                fvc::div(phi_nei,K);
                //or
                //+
                //fvc::div(phi,K);
        }
    }
    
    KSource = (p - p.oldTime())/runTime.deltaT() + fvc::div
            (
                (
                    linearInterpolate
                    (
                        (-turbulence->devRhoReff()) & U
                    ) & mesh.Sf()
                )()
            );
    


