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
    pimpleCentralDyMFoam

Description
    Pressure-based semi implicit compressible flow solver based on
    central-upwind schemes of  Kurganov and Tadmor with mesh motion and LTS
    support for steady-state calculations. Turbulence model
    is run-time selectable.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "cellQuality.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "kappaFunction.H"
#include "correctCentralACMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    pimpleControl pimple(mesh);

    #include "createRDeltaT.H"
    #include "createRDeltaTVariables.H"
    #include "createTimeControls.H"

    #include "createFields.H"
    Info << "All fields were created" << endl;
    #include "readAdditionalPimpleControl.H"
    Info << "Creating central fields" << endl;
    #include "createCommonCentralFields.H"
    Info << "Creating controls" << endl;
    #include "createCentralMeshControls.H"

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

    #include "updateCentralWeights.H"
    phi_own = phiv_own*rho_own;
    phi_nei = phiv_nei*rho_nei;
    surfaceScalarField mphi_own = phi*0.0;
    surfaceScalarField mphi_nei = mphi_own;
    #include "updateKappa.H"
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
        #include "readCentralMeshControls.H"

        {
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

            Info<< "Time = " << runTime.timeName() << nl << endl;

            // Do any mesh changes
            mesh.update();

            if (mesh.changing())
            {
                #include "updateFaceAreas.H"

                if (correctPhi)
                {
                    #include "centralCorrectPhi.H"

                    phi_own = phi_own + (1.0 - kappa) * phi_nei;
                    phi_nei = kappa * phi_nei;
                }
                else
                {
                    mphi_own = alpha_own*rho_own*mesh.phi();
                    mphi_nei = alpha_nei*rho_nei*mesh.phi();

                    //make fluxes relative
                    phi_own -= (mphi_own + (1.0 - kappa)*mphi_nei);
                    phi_nei -= (mphi_nei*kappa);
                    phi = phi_own + phi_nei;
                }
            }
        }

        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "centralMeshCourantNo.H"
            #include "markBadQualityCells.H"
        }

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
                #include "hEqn.H"
            }

            // --- Solve pressure (PISO)
            if (pimple.finalIter())
            {
                while (pimple.correct())
                {
                    #include "pressureVelocityCorrDyM.H"
                }
            }
            else
            {
                #include "pressureVelocityCorrDyM.H"
            }
            
            if (!updateEnergyInPISO)
            {
                mphi_own = alpha_own * rho_own * mesh.phi();
                mphi_nei = alpha_nei * rho_nei * mesh.phi();

                phi_own += mphi_own;
                phi_nei += mphi_nei;
                phi = phi_own + phi_nei;
                #include "updateKappa.H"
                #include "updateMechanicalFields.H"
            }

            if (!pimple.finalIter())
            {
                phi_own -= (mphi_own + (1.0 - kappa)*mphi_nei);
                phi_nei -= (mphi_nei*kappa);
                phi = phi_own + phi_nei;
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
