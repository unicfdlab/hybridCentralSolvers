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
    twoPhaseMixingCentralDyMFoam

Description
    Transient Eulerian two-phase solver with dynamic meshes. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface
    is accounted at the diffusion approximation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "kappaFunction.H"
#include "centralMULES.H"
#include "fvOptions.H"
#include "cellQuality.H"
#include "fvcSmooth.H"
#include "compressibleTwoPhaseMixtureThermo.H"
#include "dynamicFvMesh.H"
#include "correctCentralACMIInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"

    #include "createRDeltaT.H"
    #include "createRDeltaTVariables.H"
    #include "createTimeControls.H"

    #include "readTimeControls.H"
    #include "createCentralCourantNo.H"
    #include "createCentralMeshControls.H"

    #include "createCommonCentralFields.H"
    #include "createRhoHatFields.H"
    #include "createTwoPhaseFields.H"

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

    #include "markBadQualityCells.H"
    #include "readCourantType.H"
    #include "centralCompressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    surfaceScalarField mphi_own
    (
        "mphi_own",
        phi * 0.0
    );

    surfaceScalarField mphi_nei
    (
        "mphi_nei",
        mphi_own
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readAdditionalPimpleControl.H"
        #include "acousticCourantNo.H"
        #include "centralCompressibleCourantNo.H"
        #include "setDeltaT.H"

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
        // --- Move mesh and update fluxes
        {

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

                if (checkMeshCourantNo)
                {
                    #include "centralMeshCourantNo.H"
                }

                #include "markBadQualityCells.H"
             }
        }

        // --- Predict velocity for new time step
        #include "massEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            if (!updateEnergyInPISO)
            {
                #include "YLiqEqn.H"
                #include "EEqn.H"
            }

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (updateEnergyInPISO)
                {
                    #include "YLiqEqn.H"
                    #include "EEqn.H"
                }

                #include "pEqnDyM.H"

                if (updateEnergyInPISO)
                {
                    mphi_own = alpha_own*rho_own*mesh.phi();
                    mphi_nei = alpha_nei*rho_nei*mesh.phi();
                    phi_own += mphi_own;
                    phi_nei += mphi_nei;
                    phi = phi_own + phi_nei;
                    #include "updateKappa.H"
                    #include "updateMechanicalFields.H"
                }
            }


            if (!updateEnergyInPISO)
            {
                mphi_own = alpha_own*rho_own*mesh.phi();
                mphi_nei = alpha_nei*rho_nei*mesh.phi();
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
