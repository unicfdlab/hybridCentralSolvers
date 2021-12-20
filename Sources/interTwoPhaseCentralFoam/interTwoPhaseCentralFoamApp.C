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
    interTwoPhaseCentralFoam
Description
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "interTwoPhaseCentralFoam.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    bool   LTS = false;

    pimpleControl pimple(mesh);

    interTwoPhaseCentralFoam Veronika(mesh, pimple);

    #include "createTimeControls.H"
    #include "setInitialDeltaT.H"

    Info<< "\nStarting time loop\n" << endl;

    Veronika.Initialize();

    while (runTime.run())
    {

        scalarField sumPhi
        (
            fvc::surfaceSum(mag(Veronika.phiC()))().primitiveField()
        );

        CoNum = gMax((sumPhi)/mesh.V().field())*runTime.deltaTValue();

        meanCoNum =
            (gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();

        Info<< "Courant Number mean: " << meanCoNum
            << " max: " << CoNum << endl;

        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Veronika.saveOld();
        Veronika.solveRho1();
        Veronika.solveRho2();

        while (pimple.loop())
        {
            Veronika.updateK();         //Calculate values of C1_, C2_, Z1_, Z2_, K_, and phi_

            Veronika.alpha1Eqnsolve();

            Veronika.massError1();

            Veronika.massError2();

            Veronika.UEqn();          //Generate fvmatrix for UEqn (note: without grad(p_))

            Veronika.TEqnsolve();

            Veronika.Compressibility(); //Update psi1_ = molM1_/(R_ * T_) and psi2_

//            Veronika.updateK();         //Calculate values of C1_, C2_, Z1_, Z2_, K_, and phi_

            Veronika.DensityThermo();

            Veronika.speedOfSound();

            // Veronika.UpdateCentralWeights(); //Calculate fluxes (phi1_own and phi1_nei)

            // Veronika.UpdateCentralFields();  //Calculate coefficients of pEqn: phi1d_own, phi1_nei, Dp1_own, and Dp2_nei

            Veronika.UpdateCentralWeightsIndividual();

            Veronika.UpdateCentralFieldsIndividual();

            Veronika.pEqnsolve();            //Solve pEqn

            Veronika.Flux();

            Veronika.updateKappa();

            Veronika.DensityThermo();        //Update rho1_ and rho2_ through rho_ = psi_*p_

            Veronika.Density();

            Veronika.ReconstructVelocity();

            Veronika.volumeFlux();

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
