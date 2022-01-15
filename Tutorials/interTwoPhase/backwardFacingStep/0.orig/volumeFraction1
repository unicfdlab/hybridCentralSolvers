/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2012                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      volumeFraction1;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];


internalField   uniform 1;

boundaryField
{
  inlet
  {
      type            fixedValue;
      value           uniform 1;
  }
  outlet
  {
      type            zeroGradient;
  }
  walls
  {
      type            zeroGradient;
  }
  symwall
  {
      type            empty;
  }
}


// ************************************************************************* //
