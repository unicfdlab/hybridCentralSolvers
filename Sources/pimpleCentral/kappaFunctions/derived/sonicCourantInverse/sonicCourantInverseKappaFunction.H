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
    Foam::fv::sonicCourantInverseKappaFunction

Description
    Switching field kappa is calculated as 1 / (c*dt/dx)

SourceFiles
    sonicCourantInverseKappaFunction.C

\*---------------------------------------------------------------------------*/

#ifndef sonicCourantInverseKappaFunction_H
#define sonicCourantInverseKappaFunction_H

#include "kappaFunction.H"
#include "volFields.H"
#include "Switch.H"

namespace Foam
{

namespace fv
{


class sonicCourantInverseKappaFunction
: public kappaFunction
{

private:

    //- forbid default constructor
    sonicCourantInverseKappaFunction();

    //- forbid copy constructor
    sonicCourantInverseKappaFunction(const sonicCourantInverseKappaFunction& );

protected:


public:

    //-
    TypeName("sonicCourantInverse");

    //- Construct from components
    sonicCourantInverseKappaFunction
    (
	const word& name,
	const word& type,
	const dictionary& parentDict,
	const fvMesh& mesh
    );

    //-
    virtual ~sonicCourantInverseKappaFunction();

    //-
    void update();

    //-
    tmp<surfaceScalarField> kappa();

    // I/O

    //-
    virtual void writeData (Ostream& ) const;

    //-
    virtual bool read (const dictionary& dict);
};

};

};

#endif

//END-OF-FILE
