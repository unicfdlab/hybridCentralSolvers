/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
       hybridCentralSolvers | Copyright (C) 2016-2018 ISP RAS (www.unicfd.ru)
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
    rhoPimpleCentralFoam

Description
    Pressure-based semi implicit compressible viscous flow of real gas solver
    based on central-upwind schemes of Kurganov and Tadmor and LTS
    support for steady-state calculations. Turbulence model
    is run-time selectable.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "localEulerDdtScheme.H"
#include "cellQuality.H"
#include "fvOptions.H"
#include "kappaFunction.H"
#include "fvcSmooth.H"
#include "correctCentralACMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createRDeltaT.H"
    #include "createRDeltaTVariables.H"
    #include "createTimeControls.H"



    #include "createFields.H"
    #include "readAdditionalPimpleControl.H"
    #include "createCommonCentralFields.H"
    #include "createRhoHatFields.H"

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    #include "createFvOptions.H"
    #include "createMRF.H"
    #include "initContinuityErrs.H"
    #include "readCourantType.H"

    #include "markBadQualityCells.H"

    #include "psiUpdateCentralFields.H"
    #include "updateKappa.H"

    #include "updateCentralWeights.H"
    #include "createCentralCourantNo.H"

    if (!LTS)
    {
        #include "centralCompressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readAdditionalPimpleControl.H"
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "acousticCourantNo.H"
            #include "centralCompressibleCourantNo.H"
            #include "readTimeControls.H"
            #include "setDeltaT.H"
        }

        runTime++;

        psi.oldTime();
        rho.oldTime();
        p.oldTime();
        U.oldTime();
        h.oldTime();
        K.oldTime();
        rhoHat.oldTime();

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Predict density
        #include "massEqn.H"

        // --- SIMPLE-like Pressure-Velocity Coupling
        while (pimple.loop())
        {
            // --- Solve turbulence
            turbulence->correct();

            // --- Solve momentum
            #define PIMPLE_RELAX
            #include "UEqn.H"

            // --- Solve energy
            if (!updateEnergyInPISO)
            {
                #include "hEqn.H"
            }

            // --- Solve pressure (PISO)
            if (pimple.finalIter())
            {
                while (pimple.correct())
                {
                    #include "pressureVelocityCorr.H"
                }
            }
            else
            {
                #include "pressureVelocityCorr.H"
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
