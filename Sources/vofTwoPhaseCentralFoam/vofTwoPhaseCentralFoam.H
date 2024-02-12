/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.
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
Class
    Foam::singlePhaseCompressibleSystem
Description
    Compressible system with one liquid phase and one gaseous phase
SourceFiles
    singlePhaseCompressibleSystem.C
\*---------------------------------------------------------------------------*/

#ifndef interTwoPhaseCentralFoam_H
#define interTwoPhaseCentralFoam_H

#include "fvm.H"
#include "fvc.H"
#include "uniformDimensionedFields.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarField.H"
#include "kappaFunction.H"
#include "CorrectPhi.H"
#include "pimpleControl.H"
#include "physicoChemicalConstants.H"
#include "interfaceProperties.H"
#include "turbulentTransportModel.H"
#include "gravityMeshObject.H"
#include "findRefCell.H"
#include "quasiIncompressibleTwoPhaseMixture.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class interTwoPhaseCentralFoam Declaration
\*---------------------------------------------------------------------------*/

    class vofTwoPhaseCentralFoam
    :
    public IOdictionary
    {
        protected:

            pimpleControl& pimple_;

            //- Velocity
            //- MUST_READ
            volVectorField U_;

            //- Velocity defined on cell boundaries
            surfaceScalarField phi_;

            //- Pressure
            //- MUST_READ
            volScalarField p_;

            //- Piesometric pressure
            volScalarField p_rgh_;

            //- Temperature
            //- MUST_READ
            volScalarField T_;

            //- Volume fraction of the first phase
            //- MUST_READ
            volScalarField volumeFraction1_;

            //- Volume fraction of the second phase
            //- NO_READ
            volScalarField volumeFraction2_;

            //- The mixture model contains properties of
            // phases and is needed to connect to turbulence
            // model
            quasiIncompressibleTwoPhaseMixture mixture_model_;

            //-
            autoPtr<incompressible::turbulenceModel> turbulence_;

            //- the temporal change of the first phase volume fraction
            volScalarField dotVF1_;

            //- the temporal change of the second phase volume fraction
            volScalarField dotVF2_;

            //- the flux of the first phase volume fraction
            surfaceScalarField phiVF1_;

            //- the flux of the second phase volume fraction
            surfaceScalarField phiVF2_;

            //- Description of the interface
            interfaceProperties interface_;

            //- Squared speed of sound of a mixture
            volScalarField C_;

            //- An interface compressibiliy factor
            volScalarField Lambda_;

            //- The H/A matrix of the Momentum equation
            volVectorField HbyA_;

            //- The 1/A matrix of the Momentum equation
            volScalarField oneByA_;

            //- Kappa KT/PISO blending field
            surfaceScalarField kappa_;

            //- 1.0 - kappa_
            surfaceScalarField onemkappa_;

// /********************Variables for Update-Central Wieghts Function*************/

            dimensionedScalar v_zero_;

            surfaceScalarField own_;

            surfaceScalarField nei_;

            //- a volumetric flux of phase 1 in direction own
            surfaceScalarField phi1_own_;

            //- a volumetric flux of phase 1 in direction nei
            surfaceScalarField phi1_nei_;

            //- a volumetric flux of phase 2 in direction own
            surfaceScalarField phi2_own_;

            //- a volumetric flux of phase 2 in direction nei
            surfaceScalarField phi2_nei_;

            //- an interpolation of phase 1 density in direction own
            surfaceScalarField rho1_own_;

            //- an interpolation of phase 1 density in direction nei
            surfaceScalarField rho1_nei_;

            //- an interpolation of phase 2 density in direction own
            surfaceScalarField rho2_own_;

            //- an interpolation of phase 2 density in direction nei
            surfaceScalarField rho2_nei_;

// /*                             Mixture props                                  */

            //-
            surfaceScalarField alpha_own_;

            //-
            surfaceScalarField alpha_nei_;

            //-
            surfaceScalarField aSf_;

            //- mixture speed of sound interpolated to faces from owners
            surfaceScalarField Cf_own_;

            //- mixture speed of sound interpolated to faces from neighbours
            surfaceScalarField Cf_nei_;

            //-
            surfaceScalarField CfSf_own_;

            //-
            surfaceScalarField CfSf_nei_;

            //-
            surfaceScalarField amaxSf_;

// /****************************Update Central Fields*****************************/

            surfaceScalarField aphiv_own_;

            surfaceScalarField aphiv_nei_;

            surfaceScalarField phiv_own_;

            surfaceScalarField phiv_nei_;

            surfaceScalarField phiHbyA_;

            surfaceScalarField phiHbyA_own_;

            surfaceScalarField phiHbyA_nei_;

            surfaceScalarField rAUf_own_;

            surfaceScalarField rAUf_nei_;

// /*******************************Reconstruct Velocity***************************/

//             volVectorField gradp_;

// /******************************Viscosity Terms*********************************/

            volScalarField TSource1_;

            volScalarField TSource2_;

// /******************************Something Essential*****************************/

            //- Continuity error of the first phase
            volScalarField E1_;

            //- Continuity error of the first phase
            volScalarField E2_;

            //- Continuity error of a mixture
            volScalarField E_;

            surfaceScalarField vF1face_;

            surfaceScalarField vF2face_;

            volScalarField Q_;

            volScalarField dpdt_;

            volScalarField gh_;

            surfaceScalarField ghf_;

            volScalarField Wp_;

            volVectorField B_;

            surfaceScalarField phib_;

        protected:

            //- Solve for density using explicit fluxes
            void solveRho
            (
                volScalarField& rhoi,
                const surfaceScalarField& phii_own,
                const surfaceScalarField& phii_nei
            );

            //- Update mass flux components using KT scheme coefficients
            void UpdateCentralMassFluxes
            (
                const volScalarField& rhoi,
                const dimensionedScalar& rho0i,
                const volScalarField& psii,
                const volVectorField& HbyA,
                const volScalarField& rbyA,
                const surfaceScalarField& alpha_own,
                const surfaceScalarField& alpha_nei,
                const surfaceScalarField& aSf,
                surfaceScalarField& phidi_own,
                surfaceScalarField& phidi_nei,
                surfaceScalarField& phi0i_own,
                surfaceScalarField& phi0i_nei,
                surfaceScalarField& Dpi_own,
                surfaceScalarField& Dpi_nei
            );

            //- Update KT central weights using given density, velocity and
            //sound velocity
            void UpdateCentralWeights
            (
                const volScalarField& rhoi,
                const surfaceScalarField& volumeFractioni,
                const volVectorField& U,
                const volScalarField& Ci,
                surfaceScalarField& rhoi_own,
                surfaceScalarField& rhoi_nei,
                surfaceScalarField& phiv_own,
                surfaceScalarField& phiv_nei,
                surfaceScalarField& alpha_own,
                surfaceScalarField& alpha_nei,
                surfaceScalarField& aSf
            );

            //- Update value of kappa and blend all neccesary flux fields
            void kappaBlend
            (
                const surfaceScalarField& kappa,
                surfaceScalarField& flux_own,
                surfaceScalarField& flux_nei
            );

            //- Writes min and max of kappa in cells for post-processing
            void writeMaxMinKappa
            (
                const surfaceScalarField& kappa
            );

            //- Combines balance matrices from different phases
            // values of matrices m1 and m2 are multiplied
            void combineMatrices
            (
                const fvScalarMatrix& m1,
                const fvScalarMatrix& m2,
                const volScalarField& vf1,
                const volScalarField& vf2,
                fvScalarMatrix& m,
                bool removeConst = false
            );

            //- Combines balance matrices from different phases
            // values of matrices m1 and m2 are multiplied
            void combineMatrices
            (
                const fvVectorMatrix& m1,
                const fvVectorMatrix& m2,
                const volScalarField& vf1,
                const volScalarField& vf2,
                fvVectorMatrix& m,
                bool removeConst = false
            );

            wordList p_rghPatchTypes();

            void p_rghUpdatePatchFields();

        public:

            TypeName("vofTwoPhaseCentralFoam");

// * * * * * * * * * *  Constructor and Destructor Functions * * * * * * * * //
            // Constructor
            vofTwoPhaseCentralFoam(const fvMesh& mesh, pimpleControl& ctrl);

            //- Destructor
            virtual ~vofTwoPhaseCentralFoam();

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

            void saveOld();

            void CharacteristicCourant();

            scalar FlowCourant();


// * * * * * * * * * * * * * * * * Main Functions  * * * * * * * * * * * * * //

            //- Solve an equation for volumeFraction1_
            virtual void LiquidVolumeFractionSolve();

            //- Solve for density of the first phase
            void solveRho1();

            //- Solve for density of the second phase
            void solveRho2();

            //- Construct Momentum equation and extract rbyA_ and HbyA_
            void UEqn();

            //- Solve an equation for Temperature
            virtual void TEqnSolve();

            //- Solve an equation for pressure
            virtual void pEqnSolve();

            //- Estimate error of Continuity equation for the first phase
            virtual void massError1();

            //- Estimate error of Continuity equation for the second phase
            virtual void massError2();

            void TSource();

//* * * * * * * * * * * * * * * Intermidiate Functions * * * * * * * * * * * *//

            void Initialize();

//* * * * * * * * * * * * * * * * Flux Functions * * * * * * * * * * * * * *//

            void updateKappa();

// * * * * * * * * * Mixture properties and thermodynamics * * * * * * * * *//

            void DensityThermo();

            //- Updates mixture and turbulence properties
            void MixtureProperties();

            void updateLambda();

            void Compressibility();

            void speedOfSound();

//* * * * * * * * * Kurganov's coefficients Individual Phases * * * * * * * *//

            //- Performs central interpolation of various fields:
            // density, etc
            void interpolateDensities();

            void UpdateCentralWeights();

            void UpdateCentralFields();

            void CalculateMassFluxes();

    }; // End of the class


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
