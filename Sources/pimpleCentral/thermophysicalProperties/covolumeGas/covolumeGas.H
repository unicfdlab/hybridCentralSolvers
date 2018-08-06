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
    Foam::covolumeGas

Description
    Co-volume gas equation of state

SourceFiles
    covolumeGas.C
    covolumeGasI.H

\*---------------------------------------------------------------------------*/

#ifndef covolumeGas_H
#define covolumeGas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class covolumeGas;

template<class Specie>
inline covolumeGas<Specie> operator+
(
    const covolumeGas<Specie>&,
    const covolumeGas<Specie>&
);

template<class Specie>
inline covolumeGas<Specie> operator-
(
    const covolumeGas<Specie>&,
    const covolumeGas<Specie>&
);

template<class Specie>
inline covolumeGas<Specie> operator*
(
    const scalar,
    const covolumeGas<Specie>&
);

template<class Specie>
inline covolumeGas<Specie> operator==
(
    const covolumeGas<Specie>&,
    const covolumeGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const covolumeGas<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class covolumeGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class covolumeGas
:
    public Specie
{

protected:

    //- volume constant
    scalar b_;

public:

    // Constructors

        //- Construct from components
        inline covolumeGas(const Specie& sp);

        //- Construct from Istream
        covolumeGas(Istream&);

        //- Construct from dictionary
        covolumeGas(const dictionary& dict);

        //- Construct as named copy
        inline covolumeGas(const word& name, const covolumeGas&);

        //- Construct and return a clone
        inline autoPtr<covolumeGas> clone() const;

        // Selector from Istream
        inline static autoPtr<covolumeGas> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<covolumeGas> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "covolumeGas<" + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (cp - cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const covolumeGas&);
        inline void operator-=(const covolumeGas&);

        inline void operator*=(const scalar);


    // Friend operators

        friend covolumeGas operator+ <Specie>
        (
            const covolumeGas&,
            const covolumeGas&
        );

        friend covolumeGas operator- <Specie>
        (
            const covolumeGas&,
            const covolumeGas&
        );

        friend covolumeGas operator* <Specie>
        (
            const scalar s,
            const covolumeGas&
        );

        friend covolumeGas operator== <Specie>
        (
            const covolumeGas&,
            const covolumeGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const covolumeGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "covolumeGasI.H"

#ifdef NoRepository
    #include "covolumeGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
