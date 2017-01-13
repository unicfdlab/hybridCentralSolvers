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
    
    rho_pos     = fvc::interpolate(rho, pos, "reconstruct(rho)");
    rho_neg     = fvc::interpolate(rho, neg, "reconstruct(rho)");
    
    rhoU_pos    = fvc::interpolate(rho*U, pos, "reconstruct(U)");
    rhoU_neg    = fvc::interpolate(rho*U, neg, "reconstruct(U)");

    phiv_pos    = ((rhoU_pos / rho_pos) & mesh.Sf()) - meshPhi;
    phiv_neg    = ((rhoU_neg / rho_neg) & mesh.Sf()) - meshPhi;
    
    cSf_pos     = fvc::interpolate(c, pos, "reconstruct(c)")*mesh.magSf();
    cSf_neg     = fvc::interpolate(c, neg, "reconstruct(c)")*mesh.magSf();
    
    ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
    am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);
    
    a_pos       = ap/(ap - am);
    aSf         = am*a_pos;
    a_neg       = 1.0 - a_pos;
    
    aphiv_pos   = a_pos*phiv_pos - aSf;
    aphiv_neg   = a_neg*phiv_neg + aSf;
    
    phid_pos = aphiv_pos * psi_pos;
    phid_neg = aphiv_neg * psi_neg;
    
    phi_pos = rho_pos * aphiv_pos;
    phi_neg = rho_neg * aphiv_neg;
    
    Dp_pos = a_pos*rAUf;
    Dp_neg = a_neg*rAUf;
    
    for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
    
        fvScalarMatrix pcorrPosEqn
        (
            fvc::div(phi_pos)
            -
            fvm::laplacian(Dp_pos,pcorr)
        );
        
        fvScalarMatrix pcorrNegEqn
        (
            fvc::div(phi_neg)
            -
            fvm::laplacian(Dp_neg,pcorr)
        );
    
        fvScalarMatrix pcorrEqn
        (
            rho.oldTime() * (mesh.V() - mesh.V0()) / mesh.V() / runTime.deltaT()
            +
            fvm::Sp(psi/runTime.deltaT(),pcorr)
            +
            pcorrPosEqn 
            +
            pcorrNegEqn
        );
        
        pcorrEqn.solve();
        
        phi_pos += pcorrPosEqn.flux();
        phi_neg += pcorrNegEqn.flux();
    }
    
    //update pressure
    rho += psi*pcorr;
}

//
//END-OF-FILE
//

