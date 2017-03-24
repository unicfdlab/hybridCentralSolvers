{
    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );
    
    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }
    
    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );
    
    mesh.setFluxRequired(pcorr.name());
    
    tmp<volScalarField> tVr
    (
        new volScalarField
        (
            IOobject
            (
                "Vr",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimless,
            //extrapolatedCalculatedFvPatchScalarField::typeName
            zeroGradientFvPatchScalarField::typeName
        )
    );
    
    tVr.ref().primitiveFieldRef() = mesh.Vsc0()/mesh.Vsc();
    tVr.ref().correctBoundaryConditions();
    
    //
    // Is it better to recalculate aSf due to changes in mesh.Sf() ?
    //
    
    surfaceScalarField rhof = linearInterpolate(rho);
    phiv_own = (fvc::interpolate(rho*U, own, "reconstruct(U)") / rho_own & mesh.Sf()) - mesh.phi();
    phiv_nei = (fvc::interpolate(rho*U, nei, "reconstruct(U)") / rho_nei & mesh.Sf()) - mesh.phi();
    
    phi_own = rho_own * (alpha_own*phiv_own - aSf);
    phi_nei = rho_nei * (alpha_nei*phiv_nei + aSf);
    
    surfaceScalarField Dpcorr_own
    (
        "Dpcorr_own",
        alpha_own*runTime.deltaT()
    );

    surfaceScalarField Dpcorr_nei
    (
        "Dpcorr_nei",
        alpha_nei*runTime.deltaT()
    );

    //phi = phi_own + phi_nei;
    //correctUphiBCs(rho, U, phi);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn_own
        (
              fvc::div(phi_own)
            - fvm::laplacian(Dpcorr_own,pcorr)
        );
        
        fvScalarMatrix pcorrEqn_nei
        (
              fvc::div(phi_nei)
            - fvm::laplacian(Dpcorr_nei,pcorr)
        );
        
        // Solve for pcorr such that the divergence of the corrected flux
        // mathes new mesh fluxes
        fvScalarMatrix pcorrEqn
        (
            rho*(1.0 - tVr())/runTime.deltaT()
            +
            fvm::Sp(psi/runTime.deltaT(),pcorr)
            +
            pcorrEqn_own
            +
            pcorrEqn_nei
        );
        
        pcorrEqn.solve();
        
        if (pimple.finalNonOrthogonalIter())
        {
            phi_own += pcorrEqn_own.flux();
            phi_nei += pcorrEqn_nei.flux();
            phi = phi_own + phi_nei;
        }
    }
    
    if (runTime.outputTime())
    {
        pcorr.write();
    }
}

//
//END-OF-FILE
//

