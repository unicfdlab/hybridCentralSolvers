Several modifications are planned for hybrid solvers aiming at improvement
their numerical scheme:

1) the introduction of flux weighting which is applicable both for single-gas and multicomponent mixture of perfect gases;
2) the removal hard-coded coupling with gradient limiting between mass fraction equations and the energy equation;
3) the introduction of refined sonic speed interpolation procedure.

The introduction of changes is planned in next steps:

1. The modification of  reactingPimpleCentralFoam solver (done).
2. The verification and validation of reactingPimpleCentralFoam using:
 - Sod's test;
 - TwoGasShock;
 - Gases diffusion;
 - Reactive Shock Tube;
 - laminar subsonic backward facing step;
 - Ladenburg Jet?
3. Uploading of the changed version of VnV tests in www.github.com/mkraposhin/VnV repository.
4. The modification and testing of reactingPimpleCentralDyMFoam (items 1-3).
5. The modification and testing of reactingLagrangianPimpleCentralFoam (items 1-3).

Afterwards, similar changes might be introduced into pimpleCentralFoam, rhoPimpleCentralFoam, twoPhaseMixingCentralFoam and other solvers.

An additional interesting finding pertains to the application of different solvers to Sod's test: while the hybrid approach is the almost the same, different solvers show different errors (interTwoPhaseCentralFoam shows the least error even compared to rhoCentralFoam) meaning that there is an omptimal configuration that must be implmented in all solvers.

