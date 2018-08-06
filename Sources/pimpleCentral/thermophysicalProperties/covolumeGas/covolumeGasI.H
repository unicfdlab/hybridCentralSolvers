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

\*---------------------------------------------------------------------------*/

#include "covolumeGas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::covolumeGas<Specie>::covolumeGas(const Specie& sp)
:
    Specie(sp),
    b_(0.0)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::covolumeGas<Specie>::covolumeGas
(
    const word& name,
    const covolumeGas<Specie>& pg
)
:
    Specie(name, pg),
    b_(pg.b_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::covolumeGas<Specie>>
Foam::covolumeGas<Specie>::clone() const
{
    return autoPtr<covolumeGas<Specie>>(new covolumeGas<Specie>(*this));
}


template<class Specie>
inline Foam::autoPtr<Foam::covolumeGas<Specie>>
Foam::covolumeGas<Specie>::New(Istream& is)
{
    return autoPtr<covolumeGas<Specie>>(new covolumeGas<Specie>(is));
}


template<class Specie>
inline Foam::autoPtr<Foam::covolumeGas<Specie>> Foam::covolumeGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<covolumeGas<Specie>>(new covolumeGas<Specie>(dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::covolumeGas<Specie>::rho(scalar p, scalar T) const
{
//    return p/(this->R()*T);
    return 1.0 / (this->R()*T/p + b_);
}


template<class Specie>
inline Foam::scalar Foam::covolumeGas<Specie>::H(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::covolumeGas<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::covolumeGas<Specie>::S(scalar p, scalar T) const
{
    return -this->R()*log(p/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::covolumeGas<Specie>::psi(scalar p, scalar T) const
{
//    return 1.0/(this->R()*T);
    return pow(this->rho(p,T)/p,2.0) * (this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::covolumeGas<Specie>::Z(scalar p, scalar T) const
{
//    return 1;
    return p / (this->rho(p,T) * this->R() * T);
}


template<class Specie>
inline Foam::scalar Foam::covolumeGas<Specie>::CpMCv(scalar p, scalar T) const
{
    return this->R();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::covolumeGas<Specie>::operator+=(const covolumeGas<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::covolumeGas<Specie>::operator-=(const covolumeGas<Specie>& pg)
{
    Specie::operator-=(pg);
}


template<class Specie>
inline void Foam::covolumeGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::covolumeGas<Specie> Foam::operator+
(
    const covolumeGas<Specie>& pg1,
    const covolumeGas<Specie>& pg2
)
{
    return covolumeGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::covolumeGas<Specie> Foam::operator-
(
    const covolumeGas<Specie>& pg1,
    const covolumeGas<Specie>& pg2
)
{
    return covolumeGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::covolumeGas<Specie> Foam::operator*
(
    const scalar s,
    const covolumeGas<Specie>& pg
)
{
    return covolumeGas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::covolumeGas<Specie> Foam::operator==
(
    const covolumeGas<Specie>& pg1,
    const covolumeGas<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
