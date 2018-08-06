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
    Foam::fv::acceleratedSonicCourantInverseKappaFunction

Description
    Switching field kappa is calculated as 1 / (c*dt/dx)**m, where m >= 1

SourceFiles
    acceleratedSonicCourantInverseKappaFunction.C

\*---------------------------------------------------------------------------*/

#ifndef acceleratedSonicCourantInverseKappaFunction_H
#define acceleratedSonicCourantInverseKappaFunction_H

#include "kappaFunction.H"
#include "volFields.H"
#include "Switch.H"
#include "fvc.H"

namespace Foam
{

namespace fv
{


class acceleratedSonicCourantInverseKappaFunction
: public kappaFunction
{

private:

    //- forbid default constructor
    acceleratedSonicCourantInverseKappaFunction();

    //- forbid copy constructor
    acceleratedSonicCourantInverseKappaFunction(const acceleratedSonicCourantInverseKappaFunction& );

protected:

    //-
    scalar power_;

    //-
    Switch reinterpolateToFaces_;

    //-
    void reinterpolateToFaces(surfaceScalarField& kappa);

public:

    //-
    TypeName("acceleratedSonicCourantInverse");

    //- Construct from components
    acceleratedSonicCourantInverseKappaFunction
    (
	const word& name,
	const word& type,
	const dictionary& parentDict,
	const fvMesh& mesh
    );

    //-
    virtual ~acceleratedSonicCourantInverseKappaFunction();

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
