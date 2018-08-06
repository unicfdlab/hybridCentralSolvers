/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Y.Liq;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    INLET
    {
       type            uniformFixedValue;
       uniformValue    sine;
       uniformValueCoeffs
       {
           amplitude   1.0;
           scale       0.5;
           level       0.5;
           frequency   constant 33.0;
       }
//        type    fixedValue;
//        value   uniform 0;
    }
    BLADE_ROT
    {
        type            zeroGradient;
    }
    HUB_SHROUD
    {
        type            empty;
    }
    AMI_INT
    {
        type            cyclicAMI;
        value           uniform 0;
    }
    AMI_EXT
    {
        type            cyclicAMI;
        value           uniform 0;
    }
    OUTLET
    {
        type            zeroGradient;
    }
    BLADE_STAT
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
