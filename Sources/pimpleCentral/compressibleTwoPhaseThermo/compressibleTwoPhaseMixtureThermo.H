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
    
Class
    Foam::compressibleTwoPhaseMixtureThermo

Description
    A two-phase mixture thermo model that uses mass fraction or volume fraction
    when appropriates to calculate fluid mixture properties.

SourceFiles
    compressibleTwoPhaseMixtureThermo.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleTwoPhaseMixtureThermo_H
#define compressibleTwoPhaseMixtureThermo_H

#include "rhoThermo.H"
#include "simplePhasePsiThermo.H"
#include "compressibleTwoPhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class compressibleTwoPhaseMixtureThermo Declaration
\*---------------------------------------------------------------------------*/

class compressibleTwoPhaseMixtureThermo
:
    public rhoThermo,
    public compressibleTwoPhaseMixture
{
    // Private data

        //-
        dimensionedScalar pMin_;

        //- Thermo-package of each phase
        autoPtr<simplePhasePsiThermo> thermoLiq_;
        autoPtr<simplePhasePsiThermo> thermoGas_;

        //- Enthalpy/Internal energy of mixture
        volScalarField he_;

    // Private Member Functions

        //- Correct the enthalpy/internal energy field boundaries
        void heBoundaryCorrection(volScalarField& he);


public:

    //- Runtime type information
    TypeName("compressibleTwoPhaseMixtureThermo");


    // Constructors

        //- Construct from mesh
        compressibleTwoPhaseMixtureThermo
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~compressibleTwoPhaseMixtureThermo();


    // Member Functions

        //- Return thermo-package of each phase
        const simplePhasePsiThermo& thermoLiq() const
        {
            return thermoLiq_();
        }

        simplePhasePsiThermo& thermoLiq()
        {
            return thermoLiq_();
        }

        const simplePhasePsiThermo& thermoGas() const
        {
            return thermoGas_();
        }

        simplePhasePsiThermo& thermoGas()
        {
            return thermoGas_();
        }

        //- Update properties
        virtual void correct();

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const;

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const;
        
        //-
        word thermoName() const;


        // Access to thermodynamic state variables
        
            //- Molecular name
            tmp<volScalarField> W() const
            {
                return 1.0 / (YLiq() / thermoLiq().W() + YGas() / thermoGas().W());
            }

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& he()
            {
                return he_;
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const
            {
                return he_;
            }

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;


        // Fields derived from thermodynamic state variables

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity ratio []
            virtual tmp<volScalarField> CpByCpv() const;

            //- Heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpv
            (
                const scalarField& p,
                const scalarField& T,
                const label patchi
            ) const;


        // Fields derived from transport state variables

            //- Thermal conductivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappa() const;

            //- Thermal conductivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappa
            (
                const label patchi
            ) const;

            //- Effective thermal conductivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal conductivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal diffusivity of mixture [J/m/s/K]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const;

            //- Mixture density
            virtual tmp<volScalarField> rho() const;

            //-
            virtual tmp<scalarField> rho(const label patchi) const;

            //-
            const dimensionedScalar& pMin() const;

            //-
            void correctRealDensities();

            //- Total enthalpy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h();

            //- Total enthalpy [J/kg]
            virtual const volScalarField& h() const;

            //- Total enthalpy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& T,
                const label patchi
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
