
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Compressible system with one phase
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
#include "psiThermo.H"
#include "CorrectPhi.H"
#include "pimpleControl.H"
#include "physicoChemicalConstants.H"
#include "interfaceProperties.H"

#include "gravityMeshObject.H"
#include "findRefCell.H"
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

            //- Velocity estimated on cell boundaries
            surfaceScalarField phi_;

            //- Pressure
            //- MUST_READ
            volScalarField p_;

//             volScalarField p_rgh_;

            //- Temperature
            //- MUST_READ
            volScalarField T_;

            //- Specific gas constant of the first phase
            //- Set in the phaseProperties file
            dimensioned<scalar> R1_;

            //- Specific gas constant of the second phase
            //- Set in the phaseProperties file
            dimensioned<scalar> R2_;

            //- Molecular mass of the first phase
            //- Set in the phaseProperties file
            dimensioned<scalar> molM1_;

            //- Molecular mass of the second phase
            //- Set in the phaseProperties file
            dimensioned<scalar> molM2_;

            //- Compressibility coefficient of the first phase
            volScalarField psi1_;

            //- Compressibility coefficient of the second phase
            volScalarField psi2_;

            //- Density of the first phase
            volScalarField rho1_;

            //- Density of the second phase
            volScalarField rho2_;

            //- Volume fraction of the first phase
            //- MUST_READ
            volScalarField volumeFraction1_;

            //- Volume fraction of the second phase
            //- NO_READ
            volScalarField volumeFraction2_;

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

            //- Total density of a mixture
            volScalarField rho_;

            // Constant density of the first phase, i.e. rho_01
            // Set in the phaseProperties file as rho01
            dimensioned<scalar> rho01_;

            // Constant density of a phase, i.e. rho_02
            // Set in the phaseProperties file as rho02
            dimensioned<scalar> rho02_;

            //- Specific heat of the first phase
            //- Set in the phaseProperties
            dimensioned<scalar> Cp1_;

            //- Specific heat of the second phase
            //- Set in the phaseProperties
            dimensioned<scalar> Cp2_;

            //- Dynamic viscousity of the first phase
            //- Set in the phaseProperties
            dimensioned<scalar> mu1_;

            //- Dynamic viscousity of the second phase
            //- Set in the phaseProperties
            dimensioned<scalar> mu2_;

            //- Prandtl number of the first phase
            //- Set in the phaseProperties
            dimensioned< scalar > Pr1_;

            //- Prandtl number of the second phase
            //- Set in the phaseProperties
            dimensioned< scalar > Pr2_;

            //- Thermal diffusivity of the first phase
            dimensioned< scalar > alpha1_;

            //- Thermal diffusivity of the second phase
            dimensioned< scalar > alpha2_;

            //- Heat capacity ratio of the first phase
            dimensioned<scalar> gamma1_;

            //- Heat capacity ratio of the second phase
            dimensioned<scalar> gamma2_;

            //- Squared speed of sound of a mixture
            volScalarField C_;

            //- An interface compressibiliy factor
            volScalarField Lambda_;

            //- The H/A matrix of the Momentum equation
            volVectorField HbyA_;

            //- The 1/A matrix of the Momentum equation
            volScalarField oneByA_;

//             //- Volumetric flux due to body forces
//             surfaceScalarField phib_;

//             //- Kappa KT/PISO blending function
//             autoPtr<fv::kappaFunction> kappaFuncPtr_;

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

            //- mixture dynamic viscosity
            volScalarField mu_;

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

//             surfaceScalarField Dp1_own_;

//             surfaceScalarField Dp1_nei_;

//             surfaceScalarField phi2d_own_;

//             surfaceScalarField phi2d_nei_;

//             surfaceScalarField Dp2_own_;

//             surfaceScalarField Dp2_nei_;

//             surfaceScalarField phi01d_own_;

//             surfaceScalarField phi01d_nei_;

//             surfaceScalarField phi02d_own_;

//             surfaceScalarField phi02d_nei_;

//             fvScalarMatrix pEqn1_own_;

//             fvScalarMatrix pEqn1_nei_;

//             fvScalarMatrix pEqn2_own_;

//             fvScalarMatrix pEqn2_nei_;

// /*******************************Reconstruct Velocity***************************/

//             volVectorField gradp_;

// /******************************Viscosity Terms*********************************/

            fvVectorMatrix divDevRhoReff_;

//             volSymmTensorField devRhoReff1_;

//             volSymmTensorField devRhoReff2_;

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

            dimensionedScalar rho1Min_;

            dimensionedScalar rho2Min_;

            volScalarField Q_;

            volScalarField dpdt_;

//             surfaceScalarField rAUf_;

//             volScalarField gh_;

//             surfaceScalarField ghf_;

//             volScalarField Wp_;

//             surfaceScalarField Wp_own_;

//             surfaceScalarField Wp_nei_;

//             volScalarField rho0_;

//             surfaceScalarField rho0ghf_;

//             volScalarField pcorr_;

//             surfaceScalarField phi1_;

//             surfaceScalarField phi2_;

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


            // const surfaceScalarField& phiC() const
            // {
            //     return phi_;
            // }

// * * * * * * * * * * * * * * * * Main Functions  * * * * * * * * * * * * * //

            //- Solve an equation for volumeFraction1_
            virtual void alpha1Eqnsolve();

            //- Solve for density of the first phase
            void solveRho1();

            //- Solve for density of the second phase
            void solveRho2();

            //- Construct Momentum equation and extract rbyA_ and HbyA_
            void UEqn();

            //- Reconstructs velocity from pressure gradient, H and A operators
            void ReconstructVelocity();

            //- Solve an equation for Temperature
            virtual void TEqnSolve();

            //- Solve an equation for pressure
            virtual void pEqnSolve();

            //- Estimate error of Continuity equation for the first phase
            virtual void massError1();

            //- Estimate error of Continuity equation for the second phase
            virtual void massError2();

            void TSource();

            void pressureGradient();

//* * * * * * * * * * * * * * * Intermidiate Functions * * * * * * * * * * * *//

            void Initialize();

//* * * * * * * * * * * * * * * * Flux Functions * * * * * * * * * * * * * *//

            void updateKappa();

//* * * * * * * * * * * * * * * Update Dencities * * * * * * * * * * * * * * *//

            void DensityThermo();

            void Density();

//* * * * * * * * * * * * * Update Depending Variables * * * * * * * * * * *//

            virtual void updateLambda();

            void Compressibility();

            void speedOfSound();

//* * * * * * * * * Kurganov's coefficients Individual Phases * * * * * * * *//

            //- Performs central interpolation of various fields:
            // density, etc
            void interpolateDensities();

            void UpdateCentralWeights();

            void UpdateCentralFields();

            void CalculateMassFluxes();

//* * * * * * * * * * * * * * * * * Viscosity * * * * * * * * * * * * * * * *//

            void divDevRhoReff();

//* * * * * * * * * * * * * * * * * * Gravity Terms * * * * * * * * * * * * * * * *//

            void GravityCorrection();

    }; // End of the class


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
