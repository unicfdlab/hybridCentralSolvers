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
    Foam::quasiIncompressibleTwoPhaseMixture
Description
    Viscosity model for quasi incompressible two phase mixture
SourceFiles
    quasiIncompressibleTwoPhaseMixture.C
\*---------------------------------------------------------------------------*/

#ifndef quasiIncompressibleTwoPhaseMixture_H
#define quasiIncompressibleTwoPhaseMixture_H

#include "dimensionedType.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "incompressible/transportModel/transportModel.H"

namespace Foam
{

 /*---------------------------------------------------------------------------*\
              Class quasiIncompressibleTwoPhaseMixture Declaration
 \*---------------------------------------------------------------------------*/

class quasiIncompressibleTwoPhaseMixture : public transportModel
{
    //- 
    const volScalarField &volumeFraction1_;

    //-
    const volScalarField &volumeFraction2_;

    //-
    const volScalarField& p_;

    //-
    const volScalarField& T_;

    /*
     * EoS properties
     */

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

    //- Constant density of the first phase, i.e. rho_01
    // Set in the phaseProperties file as rho01
    dimensioned<scalar> rho01_;

    //- Constant density of a phase, i.e. rho_02
    // Set in the phaseProperties file as rho02
    dimensioned<scalar> rho02_;

    //- Minimum value of the first phase density
    dimensioned<scalar> rho1Min_;

    //- Minimum value of the second phase density
    dimensioned<scalar> rho2Min_;

    /*
     * Thermodynamic properties
     */

    //- Specific heat of the first phase
    //- Set in the phaseProperties
    dimensioned<scalar> Cp1_;

    //- Specific heat of the second phase
    //- Set in the phaseProperties
    dimensioned<scalar> Cp2_;

    //- Heat capacity ratio of the first phase
    dimensioned<scalar> gamma1_;

    //- Heat capacity ratio of the second phase
    dimensioned<scalar> gamma2_;

    /*
     * Transport properties
     */

    //- Dynamic viscousity of the first phase
    //- Set in the phaseProperties
    dimensioned<scalar> mu1_;

    //- Dynamic viscousity of the second phase
    //- Set in the phaseProperties
    dimensioned<scalar> mu2_;

    //- Prandtl number of the first phase
    //- Set in the phaseProperties
    dimensioned<scalar> Pr1_;

    //- Prandtl number of the second phase
    //- Set in the phaseProperties
    dimensioned<scalar> Pr2_;

    //- Thermal diffusivity of the first phase
    dimensioned<scalar> alpha1_;

    //- Thermal diffusivity of the second phase
    dimensioned<scalar> alpha2_;

    //- Turbulent Prtandtl Number
    dimensioned<scalar> Prt_;

    /*
     * Derived volume fields
     */

    //- Compressibility coefficient of the first phase
    volScalarField psi1_;

    //- Compressibility coefficient of the second phase
    volScalarField psi2_;

    //- Density of the first phase
    volScalarField rho1_;

    //- Density of the second phase
    volScalarField rho2_;

    // rho0 of the mixture
    volScalarField rho0_;

    //- Mixture density
    volScalarField rho_;

    //- Mixture dynamic viscosity
    volScalarField mu_;


public:
 
     //- Runtime type information
     TypeName("quasiIncompressibleTwoPhaseMixture");
 
 
     // Constructors
 
        //- Construct from mixture components
        // and dictionary
        quasiIncompressibleTwoPhaseMixture
        (
            const IOdictionary &dict,
            const volScalarField &volumeFraction1,
            const volScalarField &volumeFraction2,
            const volScalarField &p,
            const volScalarField &T
        );
 
 
    //- Destructor
    virtual ~quasiIncompressibleTwoPhaseMixture();
 
 
    // Member Functions
 
        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read(){return true;};

        //- Correct compressibility coefficients of phases
        void compressibility();

        //- Correct density of phases
        void density();

        //- Correct rho0 of the mixture
        void density0();

        //- Save old time values for derived fields
        void saveOldTime();

        //- limit rho1
        void limitRho1();

        //- limit rho2
        void limitRho2();

    // Access functions

    //- Returns R1
    const dimensioned<scalar> &R1() const
    {
        return R1_;
    }

    //- Returns R2
    const dimensioned<scalar> &R2() const
    {
        return R2_;
    }

    //- Returns molM1
    const dimensioned<scalar> &molM1() const
    {
        return molM1_;
    }

    //- Returns molM2
    const dimensioned<scalar> &molM2() const
    {
        return molM2_;
    }

    //- Returns rho01
    const dimensioned<scalar> &rho01() const
    {
        return rho01_;
    }

    //- Returns rho02
    const dimensioned<scalar> &rho02() const
    {
        return rho02_;
    }

    //- Returns Cp1
    const dimensioned<scalar> &Cp1() const
    {
        return Cp1_;
    }

    //- Return Cp2
    const dimensioned<scalar> &Cp2() const
    {
        return Cp2_;
    }

    //- Returns gamma1
    const dimensioned<scalar> &gamma1() const
    {
        return gamma1_;
    }

    //- Returns gamma2
    const dimensioned<scalar> &gamma2() const
    {
        return gamma2_;
    }

    //- Returns mu1
    const dimensioned<scalar> &mu1() const
    {
        return mu1_;
    }

    //- Returns mu2
    const dimensioned<scalar> &mu2() const
    {
        return mu2_;
    }

    //- Returns Pr1
    const dimensioned<scalar> &Pr1() const
    {
        return Pr1_;
    }

    //- Returns Pr2
    const dimensioned<scalar> &Pr2() const
    {
        return Pr2_;
    }

    //- Returns alpha1
    const dimensioned<scalar> &alpha1() const
    {
        return alpha1_;
    }

    //- Returns alpha2
    const dimensioned<scalar> &alpha2() const
    {
        return alpha2_;
    }

    //- Returns Prt
    const dimensioned<scalar> &Prt() const
    {
        return Prt_;
    }

    //- Returns psi1
    const volScalarField &psi1() const
    {
        return psi1_;
    }

    //- Returns psi2
    const volScalarField &psi2() const
    {
        return psi2_;
    }

    //- Returns rho1
    const volScalarField &rho1() const
    {
        return rho1_;
    }

    //- Returns rho2
    const volScalarField &rho2() const
    {
        return rho2_;
    }

    //- Returns rho1 non-const
    volScalarField &rho1()
    {
        return rho1_;
    }

    //- Returns rho2 non-const
    volScalarField &rho2()
    {
        return rho2_;
    }

    //- Returns rho0
    const volScalarField &rho0() const
    {
        return rho0_;
    }

    //- Returns rho
    const volScalarField &rho() const
    {
        return rho_;
    }

    //- Returns mu
    const volScalarField &mu() const
    {
        return mu_;
    }
};

}

#endif

