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
    pisoCentralFoam

Description
    Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

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
    
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    scalar initialDeltaT = -VGREAT;
    
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
    #include "updateKappa.H"
    #include "createCentralCourantNo.H"
    
    if (!LTS)
    {
        #include "centralCompressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
    Info << "tuSf().boundaryField().size() = " << tuSf().boundaryField().size() << endl;
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
                    
                    phi_nei += (1.0 - kappa) * phi_own;
                    phi_own *= kappa;
                }
            }
        }
        
        if (mesh.changing() && checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
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
            while (pimple.correct())
            {
                if (updateEnergyInPISO) //update each iteration before pressure
                {
                    #include "hEqn.H"
                }
                
                #include "pEqnDyM.H"
                
                if (updateEnergyInPISO)
                {
                    #define PISOCENTRALFOAM_LTS
                    
                    //// --- update weightings for central scheme
                    //#include "updateCentralWeights.H"

                    surfaceScalarField mphi_own = alpha_own * rho_own * fvc::meshPhi(rho,U);
                    surfaceScalarField mphi_nei = alpha_nei * rho_nei * fvc::meshPhi(rho,U);
                    
                    phi_nei += mphi_nei;
                    phi_own += mphi_own;
                    phi = phi_own + phi_nei;

                    // --- update blending function
                    #include "updateKappa.H"

                    phi_nei -= (mphi_nei + (1.0 - kappa)*mphi_own);
                    phi_own -= (kappa*mphi_own);
                    phi = phi_own + phi_nei;

                    // --- update mechanical fields
                    #include "updateMechanicalFields.H"
                }
            }
            

            if (!updateEnergyInPISO)
            {
                //// --- update weightings for central scheme
                //#include "updateCentralWeights.H"
                
                surfaceScalarField mphi_own = alpha_own * rho_own * fvc::meshPhi(rho,U);
                surfaceScalarField mphi_nei = alpha_nei * rho_nei * fvc::meshPhi(rho,U);
                
                phi_nei += mphi_nei;
                phi_own += mphi_own;
                phi = phi_own + phi_nei;
                
                // --- update blending function
                #include "updateKappa.H"
                
                phi_nei -= (mphi_nei + (1.0 - kappa)*mphi_own);
                phi_own -= (kappa*mphi_own);
                phi = phi_own + phi_nei;

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
