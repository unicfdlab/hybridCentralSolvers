/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    twoPhaseMixingCentralFoam

Description
    Transient Eulerian two-phase solver. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface
    is accounted at the diffusion approximation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "turbulenceModel.H"

#include "centralMULES.H"
#include "cellQuality.H"
#include "compressibleTwoPhaseMixtureThermo.H"
#include "kappaFunction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    dimensionedScalar v_zero ("v_zero", dimVolume/dimTime, 0.0);
    
    #include "createCommonCentralFields.H"
    #include "createAdditionalSurfaceFields.H"
    #include "markBadQualityCells.H"
    #include "readCourantType.H"

    #include "createCentralCourantNo.H"
    {
        #include "readTimeControls.H"
        #include "setInitialDeltaT.H"
    }
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "acousticCourantNo.H"
        #include "centralCompressibleCourantNo.H"
        #include "setDeltaT.H"
        #include "readPIMPLEControls.H"

        runTime++;
        
        thermo.he().oldTime();
        rho.oldTime();
        p.oldTime();
        psi.oldTime();
        YLiq.oldTime();
        YGas.oldTime();
        YbarLiq.oldTime();
        YbarGas.oldTime();
        rhoHat.oldTime();
        K.oldTime();

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        for (label oCorr = 0; oCorr < nOuterCorr; oCorr++)
        {
            #include "MixtureRhoEqn.H"
            
            scalarField allFacesLambda(mesh.nFaces(), 1.0);
            slicedSurfaceScalarField lambdaCoeffs
            (
                IOobject
                (
                    "lambdaCoeffs",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimless,
                allFacesLambda,
                false   // Use slices for the couples
            );

            #include "YLiqEqn.H"
            #include "UEqn.H"
            #include "EEqn.H"

            // --- Pressure corrector loop
            for (label corr = 1; corr <= nCorr; corr++)
            {
                #include "pEqn.H"
            }
            
            #include "updateKappa.H"
            
            dpdt = fvc::ddt(p);
            KChange = fvc::ddt(rho,K) + fvc::div(phi_pos,K) + fvc::div(phi_neg,K)
                - fvc::div( ((-turbulence->devRhoReff()) & U) );
            
            turbulence->correct();
        }

        if(runTime.write())
        {
            c.write();
            YbarLiq.write();
        }
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
