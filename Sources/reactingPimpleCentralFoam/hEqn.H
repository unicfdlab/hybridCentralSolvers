{
    fvScalarMatrix hConvImpl
    (
        fvm::ddt(rho,h)
        +
        fvm::div(phi_own,h)
        +
        fvm::div(phi_nei,h)
    );

    /*
     *
     * Diffusive fluxes of energy due to mass diffusion
     *
     */
    PtrList<volScalarField> hi (Y.size());
    forAll (hi, iCmpt)
    {
        hi.set
        (
            iCmpt,
            new volScalarField
            (
                IOobject
                (
                    "h" + Y[iCmpt].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimEnergy/dimMass
            )
        );
        
        scalarField& hiIF = hi[iCmpt].ref();
        const scalarField& pIF  = p.internalField();
        const scalarField& TIF  = T.internalField();
        
        forAll(hiIF, iCell)
        {
            hiIF[iCell] = thermo.composition().Hs(iCmpt, pIF[iCell], TIF[iCell]);
        }
    
        forAll(hi[iCmpt].boundaryField(), iPatch)
        {
            fvPatchScalarField& hip = hi[iCmpt].boundaryFieldRef()[iPatch];
            const fvPatchScalarField& pp = p.boundaryField()[iPatch];
            const fvPatchScalarField& Tp = T.boundaryField()[iPatch];
            forAll(hip, iFace)
            {
                hip[iFace] = thermo.composition().Hs(iCmpt, pp[iFace], Tp[iFace]);
            }
        }
    }

    volScalarField alphahEff(turbulence->alphaEff());
    surfaceScalarField dzetaPhi
    (
        "dzetaPhi",
        //contribution from relative motion of components (diffusion in Yi transport)
        fvc::flux
        (
            diffusiveFlux[0],
            hi[0],
            "div(rhoi*Uri,hi)"
        )
        +
        // contribution due to casting diffusion term from temperature to system enthalpy
        fvc::flux
        (
            fvc::snGrad(Y[0]) * mesh.magSf(),
            hi[0] * alphahEff,
            "div(alphaEff*hi,grad(Yi))"
        )
    );
    
    for (label iCmpt = 1; iCmpt < Y.size(); iCmpt++)
    {
        dzetaPhi +=
            //contribution from relative motion of components (diffusion in Yi transport)
            fvc::flux
            (
                diffusiveFlux[iCmpt],
                hi[iCmpt],
                "div(rhoi*Uri,hi)"
            )
            +
            // contribution due to casting diffusion term from temperature to system enthalpy
            fvc::flux
            (
                fvc::snGrad(Y[iCmpt]) * mesh.magSf(),
                hi[iCmpt] * alphahEff,
                "div(alphaEff*hi,grad(Yi))"
            );
    }

    #include "hExtSources.H"

    fvScalarMatrix hEqn
    (
        hConvImpl + fvc::div(dzetaPhi)
        +
        KChange
        ==
        fvm::laplacian(turbulence->alphaEff(), h)
        +
        KSource
        +
        hExtSources
        +
        fvOptions(rho,h)
    );
    
    fvOptions.constrain(hEqn);
    
    hEqn.solve();
    
    fvOptions.correct(h);
    
    thermo.correct();
    
    Info << "max/min T: " << max(T).value() << "/" << min(T).value() << endl;
}

