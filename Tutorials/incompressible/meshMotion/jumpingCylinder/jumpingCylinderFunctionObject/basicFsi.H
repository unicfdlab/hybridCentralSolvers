/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::basicFsi

Group
    grpbasicFsiFunctionObjects

Description
    

SeeAlso
    

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef basicFsi_H
#define basicFsi_H


#include "forces.H"
#include "volFieldsFwd.H"
#include "Tuple2.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class basicFsi Declaration
\*---------------------------------------------------------------------------*/

class basicFsi
:
    public forces
{

protected:

        //- cylinder mass
        scalar M_;
        
        //- damping coefficient
        scalar C_;
        
        //- rigidity coefficient
        scalar K_;
        
        //- ratio of real cylinder length to depth of domain
        scalar R_;
        
        //- maximum amplitude of displacement
        scalar Ymax_;
        
        //- current state of system (y,Vy)
        Pair<scalar> Y_; 
        
        //- old state of system (y, Vy)
        Pair<scalar> Yold_;
        
        //- Pointer to file with results
        autoPtr<OFstream> fsiResFilePtr_;

protected:

        //- Disallow default bitwise copy construct
        basicFsi(const basicFsi&);

        //- Disallow default bitwise assignment
        void operator=(const basicFsi&);
        
        //-
        void createFsiOutFile(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("basicFsi");


    // Constructors

        //- Construct from given Time and dictionary.
        basicFsi
        (
            const word& name,
            const Time&,
            const dictionary&
        );
        
        //- Construct from given objectRegistry and dictionary.
        basicFsi
        (
            const word& name,
            const objectRegistry&,
            const dictionary&
        );

    //- Destructor
    virtual ~basicFsi();


    // Member Functions

        //- Read the basicFsi data
        virtual bool read(const dictionary&);
        
        //-
        virtual bool execute();
        
        //- Write the basicFsi
        virtual bool write();
        
        //- Distributes displacements between slave processes
        // and sets cellDisplacement field Y component on patch
        void setDisplacements(volVectorField& yDispl);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
