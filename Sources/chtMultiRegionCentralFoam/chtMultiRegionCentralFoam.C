/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    chtMultiRegionCentralFoam

Description
    Pressure-based semi implicit solver, based on hybrid central-upwind schemes
    of Kurganov and Tadmor for conjugate simulation of compressible flows of 
    perfect gas (Mach number is ranging from 0 to 6) and solid body heat
    transfer.
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "localEulerDdtScheme.H"
#include "cellQuality.H"
#include "fvOptions.H"
#include "kappaFunction.H"
#include "fvcSmooth.H"
#include "correctCentralACMIInterpolation.H"
#include "regionProperties.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "coordinateSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMeshes.H"
    
    pimpleControl pimple(fluidMesh);
    
    #include "createRDeltaT.H"
    #include "createRDeltaTVariables.H"
    #include "createTimeControls.H"
    
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
        #include "readAdditionalPimpleControl.H"
        #include "acousticCourantNo.H"
        #include "centralCompressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"
        
        runTime++;
        
        psi.oldTime();
        rho.oldTime();
        p.oldTime();
        U.oldTime();
        h.oldTime();
        K.oldTime();
        
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Predict density
        #include "massEqn.H"
        
        // --- SIMPLE-like Pressure-Velocity Coupling
        while (pimple.loop())
        {
            // --- Solve turbulence
            turbulence->correct();
            
            // --- Solve momentum
            #include "UEqn.H"
            
            // --- Solve energy
            if (!updateEnergyInPISO)
            {
                #include "solveSolidRegions.H"
                #include "hEqn.H"
            }
            
            // --- Solve pressure (PISO)
            while (pimple.correct())
            {
                if (updateEnergyInPISO) //update each iteration before pressure
                {
                    #include "solveSolidRegions.H"
                    #include "hEqn.H"
                }
                
                #include "pEqn.H"
                if (updateEnergyInPISO)
                {
                    // --- update blending function
                    #include "updateKappa.H"
                    
                    // --- update mechanical fields
                    #include "updateMechanicalFields.H"
                }
            }
            
            if (!updateEnergyInPISO)
            {
                // --- update blending function
                #include "updateKappa.H"
                
                // --- update mechanical fields
                #include "updateMechanicalFields.H"
            }
        }

        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
