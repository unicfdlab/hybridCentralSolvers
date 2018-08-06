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
    Foam::simplePhasePsiThermo

Description
    A simple thermo model for phase of equilibrium multiphase mixture

SourceFiles
    simplePhasePsiThermo.C

\*---------------------------------------------------------------------------*/

#ifndef simplePhasePsiThermo_H
#define simplePhasePsiThermo_H

#include "fvMesh.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"

namespace Foam
{

class simplePhasePsiThermo
{

public:

    scalar Cp(const scalar& p, const scalar& T)
    {
	return Cp_;
    }

private:

    //-
    simplePhasePsiThermo();

    //-
    simplePhasePsiThermo(const simplePhasePsiThermo&);

protected:

    //-
    word name_;

    //-
    const fvMesh& mesh_;

    //-
    const volScalarField& p_;

    //-
    const volScalarField& T_;

    //-
    volScalarField he_;

    //-
    volScalarField psi_;

    //-
    volScalarField alpha_;

    //-
    volScalarField rho_;

    //-
    volScalarField rhoZero_;

    //thermophysical and transport data

    //- R/W
    scalar R_;

    //- isobaric heat capacity
    scalar Cp_;

    //- isochoric heat capacity
    scalar Cv_;

    //- dynamic viscosity
    scalar mu_;

    //- Pr
    scalar Pr_;

    //-
    scalar rho0_;

    //-
    scalar p0_;

protected:

    //-
    wordList heBoundaryBaseTypes();

    //-
    wordList heBoundaryTypes();

    //-
    void heBoundaryCorrection(volScalarField& h);

public:

    //-
    simplePhasePsiThermo(const fvMesh& mesh, const dictionary& dict);

    //-
    virtual ~simplePhasePsiThermo();

    //-
    virtual void correct();

    // Access to thermodynamic state variables

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

	//-
	dimensionedScalar mu() const
	{
	    return dimensionedScalar (name_ + "mu", dimMass / dimTime / dimLength, mu_);
	}
	
	//-
	dimensionedScalar W() const
	{
	    return Foam::constant::physicoChemical::R / this->R_;
	}

	//-
	const volScalarField& psi() const
	{
	    return psi_;
	}

	//-
	const volScalarField& alpha() const
	{
	    return alpha_;
	}

	volScalarField& rho()
	{
	    return rho_;
	}

	volScalarField& rhoZero()
	{
	    return rhoZero_;
	}

        void correctDensity(const volScalarField& pLimited)
        {
            forAll(rho_, celli)
            {
                rho_[celli] = rhoZero_[celli] + psi_[celli] * pLimited[celli];
            }

            forAll(T_.boundaryField(), patchi)
            {
                fvPatchScalarField& prho = rho_.boundaryFieldRef()[patchi];
                const fvPatchScalarField& ppsi = psi_.boundaryField()[patchi];
                const fvPatchScalarField& pp = pLimited.boundaryField()[patchi];

                prho = rhoZero_.boundaryField()[patchi] +
                    ppsi*pp;
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

//
//END-OF-FILE
//
