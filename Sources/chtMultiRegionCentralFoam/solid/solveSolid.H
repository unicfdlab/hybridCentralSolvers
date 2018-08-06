if (pimple.finalIter())
{
    solidMesh.data::add("finalIteration", true);
}

{
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(solidBetav,solidRho,solidH)
          - (
               solidThermo.isotropic()
             ? fvm::laplacian(solidBetav*solidThermo.alpha(), solidH, "laplacian(alpha,h)")
             : fvm::laplacian(solidBetav*taniAlpha(), solidH, "laplacian(alpha,h)")
            )
          ==
            solidFvOptions(solidRho, solidH)
        );

        hEqn.relax();

        solidFvOptions.constrain(hEqn);

        hEqn.solve(solidMesh.solver(solidH.name()));

        solidFvOptions.correct(solidH);
    }
}

solidThermo.correct();

Info<< "Min/max T:" << min(solidThermo.T()).value() << ' '
    << max(solidThermo.T()).value() << endl;

if (pimple.finalIter())
{
    solidMesh.data::remove("finalIteration");
}
