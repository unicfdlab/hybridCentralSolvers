forAll(solidRegions, i)
{
    Info<< "\nSolving for solid region "
        << solidRegions[i].name() << endl;
    #include "setRegionSolidFields.H"
    #include "solveSolid.H"
}
