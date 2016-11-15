MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
    fvm::ddt(rho,U)
    +
    fvm::div(phi_own,U)
    +
    fvm::div(phi_nei,U)
    +
    MRF.DDt(rho,U)
    +
    turbulence->divDevRhoReff(U)
    ==
    fvOptions(rho, U)
);

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
	UEqn == -gradp
    );
    
    fvOptions.correct(U);
}
