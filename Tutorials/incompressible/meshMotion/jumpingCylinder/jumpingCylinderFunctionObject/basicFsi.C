/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "basicFsi.H"
#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(basicFsi, 0);
    
    addToRunTimeSelectionTable(functionObject, basicFsi, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::basicFsi::basicFsi
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces
    (
        name,
        runTime,
        dict
    ),
    M_(0.0),
    C_(0.0),
    K_(0.0),
    R_(0.0),
    Ymax_(0.0),
    Y_ (0.0, 0.0),
    Yold_(0.0, 0.0),
    fsiResFilePtr_(0)
{
    this->read(dict);
    this->createFsiOutFile(dict);
}

Foam::functionObjects::basicFsi::basicFsi
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forces
    (
        name,
        obr,
        dict
    ),
    M_(0.0),
    C_(0.0),
    K_(0.0),
    R_(0.0),
    Ymax_(0.0),
    Y_ (0.0, 0.0),
    Yold_(0.0, 0.0),
    fsiResFilePtr_(0)
{
    this->read(dict);
    this->createFsiOutFile(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::basicFsi::~basicFsi()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::basicFsi::createFsiOutFile(const dictionary& dict)
{
    if (Pstream::master())
    {
        fsiResFilePtr_.reset
        (
            new OFstream(dict.lookup("results"))
        );
        
        fsiResFilePtr_() << "Time;Y;Vy;Fy" << endl;
    }
}

bool Foam::functionObjects::basicFsi::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }
    
    dict.lookup("M") >> M_;
    
    dict.lookup("C") >> C_;
    
    dict.lookup("K") >> K_;
    
    dict.lookup("R") >> R_;
    
    dict.lookup("Ymax") >> Ymax_;
    
    return true;
}
void Foam::functionObjects::basicFsi::setDisplacements(volVectorField& yDispl)
{
    if (Pstream::parRun())
    {
        Pstream::scatter<scalar>(Y_.first());
    }
    
    vector YPatch (0.0, Y_.first(), 0.0);
    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
	label patchId = iter.key();
	forAll(yDispl.boundaryField()[patchId], faceI)
	{
	    yDispl.boundaryFieldRef()[patchId][faceI] = YPatch;
	}
    }
}

bool Foam::functionObjects::basicFsi::execute()
{
    return true;
}

bool Foam::functionObjects::basicFsi::write()
{

    if (!forces::execute())
    {
        return false;
    }

    volVectorField& yDispl = 
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>("cellDisplacement")
        );

    if (Pstream::master())
    {
        scalar dt = yDispl.mesh().time().deltaT().value();
        scalar ct = yDispl.mesh().time().value();
        
        vector force = forceEff();
        scalar yForce = force.y();
        
        //Runge-Kutta 2-nd order method
        Pair<scalar> Ymid;
        
        Ymid.first() = Yold_.first() + 0.5*dt*Yold_.second();
        Ymid.second()= Yold_.second() + 0.5*dt*(-C_*Yold_.second() - K_*Yold_.first() + R_*yForce) / M_;
        
        Y_.first() = Yold_.first() + dt*Ymid.second();
        Y_.second()= Yold_.second() + dt*(-C_*Ymid.second() - K_*Ymid.first() + R_*yForce) / M_;
        
        if (mag(Y_.first()) >= Ymax_)
        {
            Y_.first() = sign(Y_.first())*Ymax_;
            Y_.second() = (Y_.first() - Yold_.first()) / dt;
        }
        
        Yold_ = Y_;
        
        Log << "yForce = " << yForce << endl;
        Log << "Y= " << Y_.first() << endl;
        Log << "Vy= " << Y_.second() << endl;
        
        fsiResFilePtr_() << ct << ";" << Y_.first() << ";" << Y_.second() << ";" << yForce << endl;

    }
    
    setDisplacements(yDispl);
    
    return true;
}



// ************************************************************************* //
