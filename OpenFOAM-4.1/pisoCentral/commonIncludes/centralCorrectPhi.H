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
    phi_own = fvc::interpolate(rho, own, "reconstruct(rho)") * 
            (
                alpha_own * 
                (
                    ((fvc::interpolate(tVr() * U * rho, own, "reconstruct(U)") / rho_own) & mesh.Sf()) - 
                    fvc::meshPhi(rho,U)
                )
                - aSf
            );

    phi_nei = fvc::interpolate(rho, nei, "reconstruct(rho)") *
            (
                alpha_nei * 
                (
                    ((fvc::interpolate(tVr() * U * rho, nei, "reconstruct(U)") / rho_nei) & mesh.Sf()) - 
                    fvc::meshPhi(rho,U)
                ) 
                + aSf
            );
    
    surfaceScalarField mphid_own
    (
        "mphid_own",
        psi_own *
        (
            alpha_own * fvc::meshPhi(rho,U)
            + aSf
        )
    );

    surfaceScalarField mphid_nei
    (
        "mphid_nei",
        psi_nei *
        (
            alpha_nei * fvc::meshPhi(rho,U)
            - aSf
        )
    );
    
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
            - fvm::div(mphid_own,pcorr) - fvm::laplacian(Dpcorr_own,pcorr)
            + fvc::div(phi_own)
        );
        
        fvScalarMatrix pcorrEqn_nei
        (
            - fvm::div(mphid_nei,pcorr) - fvm::laplacian(Dpcorr_nei,pcorr)
            + fvc::div(phi_nei)
        );
        
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divRhoU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            rho*(1.0 - tVr)/runTime.deltaT()
            +
            fvm::Sp(psi/runTime.deltaT(),pcorr)
            //+
            //fvc::div(phi)
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

