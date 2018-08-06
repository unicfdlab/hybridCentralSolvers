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
    Foam::compressibleTwoPhaseMixture

Description
    A two-phase mixture model with both mass fraction and volume fraction
    and a const reference to the effective density.

SourceFiles
    compressibleTwoPhaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleTwoPhaseMixture_H
#define compressibleTwoPhaseMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class compressibleTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class compressibleTwoPhaseMixture
{
protected:

    // Protected data

        //-
        word liqPhaseName_;

        //-
        word gasPhaseName_;

        //-
        volScalarField YLiq_;

        //-
        volScalarField YGas_;

        //-
        volScalarField YbarLiq_;

        //-
        volScalarField YbarGas_;

        //-
        volScalarField rhoEff_;


public:

    // Constructors

        //- Construct from components
        compressibleTwoPhaseMixture
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~compressibleTwoPhaseMixture()
    {}


    // Member Functions

        //- Update volume fractions
        void updateVolFrac(const volScalarField& rhoLiq, const volScalarField& rhoGas);

        const word& liqPhaseName() const
        {
            return liqPhaseName_;
        }

        const word& gasPhaseName() const
        {
            return gasPhaseName_;
        }

        //- Return mass fraction of each phase
        const volScalarField& YLiq() const
        {
            return YLiq_;
        }

        volScalarField& YLiq()
        {
            return YLiq_;
        }

        const volScalarField& YGas() const
        {
            return YGas_;
        }

        volScalarField& YGas()
        {
            return YGas_;
        }

        //- Return volume fraction of each phase
        const volScalarField& YbarLiq() const
        {
            return YbarLiq_;
        }

        volScalarField& YbarLiq()
        {
            return YbarLiq_;
        }

        const volScalarField& YbarGas() const
        {
            return YbarGas_;
        }

        volScalarField& YbarGas()
        {
            return YbarGas_;
        }

        //- Return const reference to effective density
        const volScalarField& rhoEff() const
        {
            return rhoEff_;
        }

        volScalarField& rhoEff()
        {
            return rhoEff_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
