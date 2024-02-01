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
#include "turbulentFluidThermoModel.H"
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
            Veronika.updateLambda();         //Calculate values of C1_, C2_, Z1_, Z2_, Lambda_, and phi_

            Veronika.alpha1Eqnsolve();
            
            Veronika.Density();

            Veronika.UEqn();             //Generate fvmatrix for UEqn (note: without grad(p_))

            Veronika.TEqnSolve();

            Veronika.Compressibility(); //Update psi1_ = molM1_/(R_ * T_) and psi2_

            Veronika.DensityThermo();

            Veronika.speedOfSound();

            Veronika.UpdateCentralWeights();

            Veronika.UpdateCentralFields();

            // Solve pressure equation
            Veronika.pEqnSolve();

            // Veronika.Flux();

            // Update rho1_ and rho2_ through rhoi_ = psii_*p_ + \rhoi_0_
            Veronika.DensityThermo();

            Veronika.Density();

            Veronika.interpolateDensities();

            Veronika.CalculateMassFluxes();

            Veronika.massError1();

            Veronika.massError2();

            Veronika.updateKappa();

            // Veronika.ReconstructVelocity();

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
