{
    dimensionedScalar rAUf ("rAUf", dimTime, 1);
    
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
    
    phi = fvc::interpolate(rho*U) & mesh.Sf();
    fvc::makeRelative(phi, rho, U);
    
    for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            rho.oldTime() * (mesh.V() - mesh.V0()) / mesh.V() / runTime.deltaT()
            +
            fvm::Sp(psi/runTime.deltaT(),pcorr)
            +
            fvc::div(phi)
            -
            fvm::laplacian(rAUf,pcorr)
        );
        
        pcorrEqn.solve();
        
        phi += pcorrEqn.flux();
    }
    
    //update pressure
    rho += psi*pcorr;
}

//
//END-OF-FILE
//

