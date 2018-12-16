    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    
    surfaceScalarField own
    (
        IOobject
        (
            "own",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("own", dimless, 1.0)
    );
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("nei", dimless, -1.0)
    );
    
    tmp<surfaceScalarField> tuMagSf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "uMagSf",
                runTime.timeName(),
                mesh
            ),
            mesh.magSf() * 1.0
        )
    );
    
    tmp<surfaceVectorField> tuSf
    (
        new surfaceVectorField
        (
            IOobject
            (
                "uSf",
                runTime.timeName(),
                mesh
            ),
            mesh.Sf() * 1.0
        )
    );
    correctCentralACMIInterpolation(tuMagSf.ref());
    correctCentralACMIInterpolation(tuSf.ref());
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*U) & tuSf()
    );
    
    surfaceScalarField rho_own
    (
        "rho_own",
        fvc::interpolate(rho, own, "reconstruct(rho)")
    );
    surfaceScalarField rho_nei
    (
        "rho_nei",
        fvc::interpolate(rho, nei, "reconstruct(rho)")
    );
    
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(rho_nei);

    surfaceScalarField psi_own
    (
        "psi_own",
        fvc::interpolate(psi, own, "reconstruct(psi)")
    );
    surfaceScalarField psi_nei
    (
        "psi_nei",
        fvc::interpolate(psi, nei, "reconstruct(psi)")
    );
    surfaceScalarField psif
    (
        "psif",
        linearInterpolate(psi)
    );

    surfaceVectorField psiU_own
    (
        "rhoU_own",
        fvc::interpolate(psi*U, own, "reconstruct(U)")
    );
    surfaceVectorField psiU_nei
    (
        "rhoU_nei",
        fvc::interpolate(psi*U, nei, "reconstruct(U)")
    );

    surfaceScalarField p_own
    (
        "p_own",
        fvc::interpolate(p, own, "reconstruct(p)")
    );
    surfaceScalarField p_nei
    (
        "p_nei",
        fvc::interpolate(p, nei, "reconstruct(p)")
    );
    
    //Update fields at ACMI interfaces
    correctCentralACMIInterpolation(psi_nei);
    
    surfaceScalarField phiv_own
    (
        "phiv_own",
        (psiU_own / psi_own)  & mesh.Sf()
    );
    surfaceScalarField phiv_nei
    (
        "phiv_nei",
        (psiU_nei / psi_nei) & mesh.Sf()
    );
    
    surfaceScalarField cf_own
    (
        "cf_own",
        fvc::interpolate(c, own, "reconstruct(psi)")
    );
    surfaceScalarField cf_nei
    (
        "cf_nei",
        fvc::interpolate(c, nei, "reconstruct(psi)")
    );
    correctCentralACMIInterpolation(cf_nei);

    surfaceScalarField cSf_own
    (
        "cSf_own",
        cf_own * mesh.magSf()
    );
    surfaceScalarField cSf_nei
    (
        "cSf_nei",
        cf_nei * mesh.magSf()
    );
    cSf_own.setOriented(true);
    cSf_nei.setOriented(true);
    
    surfaceScalarField ap
    (
        "ap",
        max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero)
    );
    surfaceScalarField am
    (
        "am",
        min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero)
    );

    surfaceScalarField alpha_own
    (
        "alpha_own",
        ap/(ap - am)
    );
    surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
    amaxSf.setOriented(true);
    
    surfaceScalarField aSf
    (
        "aSf",
        am*alpha_own
    );
    surfaceScalarField alpha_nei
    (
        "alpha_nei",
        1.0 - alpha_own
    );

    if (fluxScheme == "Kurganov")
    {
        //do nothing, field is updated at each time step
    }
    else if (fluxScheme == "Tadmor")
    {
        //fields are equal
        alpha_own = 0.5;
        alpha_nei = 1.0 - alpha_own;
    }
    else
    {
        
    }
    
    phiv_own *= alpha_own;
    phiv_nei *= alpha_nei;
    
    surfaceScalarField aphiv_own("aphiv_own", phiv_own - aSf);
    surfaceScalarField aphiv_nei("aphiv_nei", phiv_nei + aSf);
    
    surfaceScalarField phid_own
    (
        "phid_own",
        aphiv_own * psi_own
    );
    surfaceScalarField phid_nei
    (
        "phid_nei",
        aphiv_nei * psi_nei
    );

    surfaceScalarField phi_own
    (
        "phi_own",
        phi * 1.0
    );
    surfaceScalarField phi_nei
    (
        "phi_nei",
        phi * 0.0
    );

    surfaceScalarField Dp_own
    (
        "Dp_own",
        alpha_own * dimensionedScalar ("zeroA", dimTime, 0.0)
    );
    surfaceScalarField Dp_nei
    (
        "Dp_nei",
        Dp_own
    );
    
    /*
     *
     * Flux scheme switching surface field
     *
     */
    autoPtr<fv::kappaFunction> kappaFuncPtr
    (
        fv::kappaFunction::New
        (
            "kappaFunction",
            mesh.solutionDict().subDict("PIMPLE").subDict("kappaFunction"),
            mesh
        )
    );
    surfaceScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    );
    forAll(kappa, iFace)
    {
        kappa[iFace] = 0.0;
    }
    forAll(kappa.boundaryField(), iPatch)
    {
        kappa.boundaryFieldRef()[iPatch] = 0.0;
    }

