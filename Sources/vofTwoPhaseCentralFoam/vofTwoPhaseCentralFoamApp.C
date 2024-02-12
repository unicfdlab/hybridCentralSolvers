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
    vofTwoPhaseCentralFoam
Description
    Solves for motion of 2 compressible (one gas and one liquid) phases
    using ACID and hybrid KNP/PIMPLE numerical scheme
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "vofTwoPhaseCentralFoam.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    bool   LTS = false;

    pimpleControl pimple(mesh);

    vofTwoPhaseCentralFoam Veronika(mesh, pimple);

    #include "createTimeControls.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    Veronika.Initialize();

    while (runTime.run())
    {

        CoNum = Veronika.FlowCourant();

        Info<< "Courant Number max: " << CoNum
            << endl;

        #include "readTimeControls.H"
        #include "setDeltaT.H"
                
        Veronika.CharacteristicCourant();

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Veronika.saveOld();
        Veronika.solveRho1();
        Veronika.solveRho2();

        while (pimple.loop())
        {
            /*
             * Calculate compressibility at the interface
             * for the advection of liquid (volumeFraction1)
             * volume fraction field
             * Calculates values of: C1_, C2_, Z1_, Z2_, Lambda_
             */
            Veronika.updateLambda();

            /*
             * Solve for evolution (advection) of liquid(volumeFraction1)
             * volume fraction field
             * Calculates:
             * - volumeFraction1_
             * - volumeFraction2_
             * - fluxes of volumeFraction1_ volumeFraction2_
             * - temporal changes of phases volume fractions
             * - face interpolations of phases volume fractions
             */
            Veronika.LiquidVolumeFractionSolve();
            

            /*
             * Update mixture density and viscosity using new composition
             * of the mixture
             * Updates turbulence properties
             */
            Veronika.MixtureProperties();

            /*
             * Create discretized equation for the mixture velocity
             */
            Veronika.UEqn();

            /*
             * Solve energy equation for the mixture temperature
             */
            Veronika.TEqnSolve();

            /*
             * Calculate compressibilities of the mixture phases
             * psi1_ and psi2_
             */
            Veronika.Compressibility();

            /*
             * Update densities of the mixture phases
             */
            Veronika.DensityThermo();

            /*
             * Update approximation of the mixture speed of sound
             */
            Veronika.speedOfSound();

            /*
             * Update weights of KNP central scheme for hyperbolic
             * parts of the system
             */
            Veronika.UpdateCentralWeights();
            Veronika.UpdateCentralFields();

            /*
             * Solve pressure equation
             */
            Veronika.pEqnSolve();

            /*
             * Update phases densities according to the new pressure field
             */
            Veronika.DensityThermo();

            /*
             * Interpolate phases densities to faces
             */
            Veronika.interpolateDensities();

            /*
             * Calculate mass fluxes from interpolated densities
             * and volumetric fluxes obtained from pressure equation
             */
            Veronika.CalculateMassFluxes();

            /*
             * Calculate mass error for the phase 1
             */
            Veronika.massError1();

            /*
             * Calculate mass error for the phase 2
             */
            Veronika.massError2();

            /*
             * Calculate blending KNP/PIMPLE blending coefficient.
             * Apply the blending coefficient to fluxes (volumetric
             * and mass).
             */
            Veronika.updateKappa();

            /*
             * Calculate mechanical energy transport terms
             */
            Veronika.TSource();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
