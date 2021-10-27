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
    Foam::fv::kappaFunction

Description
    Abstract base class for the estimation of kappa surface field
    which is used to switch between Kurganov-Tadmor convective fluxes
    and standard PIMPLE fluxes

SourceFiles
    kappaFunction.C

\*---------------------------------------------------------------------------*/
#ifndef kappaFunction_H
#define kappaFunction_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "HashTable.H"
#include "Switch.H"


namespace Foam
{

class fvMesh;
class Time;

namespace fv
{


class kappaFunction
{

private:

    //- forbid default constructor
    kappaFunction();

    //- forbid copy constructor
    kappaFunction(const kappaFunction& );

protected:

    //-
    const word name_;

    //-
    const word scType_;

    //-
    const fvMesh& mesh_;

    //-
    const Time& runTime_;

    //-
    const dictionary coeffs_;

    //-
    Switch writeMaxMin_;

    //-
    Switch resetCoupledBoundaries_;

protected:

    //-
    void writeMaxMinKappa (const surfaceScalarField& kappa);

    //-
    void resetCoupledBoundaries(surfaceScalarField& kappa);

public:

    //-
    TypeName("kappaFunction");

    declareRunTimeSelectionTable
    (
        autoPtr,
        kappaFunction,
        dictionary,
        (
            const word& name,
            const word& type,
            const dictionary& parentDict,
            const fvMesh& mesh
        ),
        (name, type, parentDict, mesh)
    );

    //- Construct from components
    kappaFunction
    (
        const word& name,
        const word& type,
        const dictionary& parentDict,
        const fvMesh& mesh
    );

    //- Return clone
    autoPtr<kappaFunction> clone() const
    {
        notImplemented("autoPtr<kappaFunction> clone() const");
        return autoPtr<kappaFunction>(nullptr);
    }

    //- Return pointer to new kappaFunction object created
    //  on the freestore from an Istream
    class iNew
    {

        //- Reference to the mesh database
        const fvMesh& mesh_;
        const word& name_;

    public:
        iNew
        (
            const fvMesh& mesh,
            const word& name
        )
        :
        mesh_(mesh),
        name_(name)
        {}

        autoPtr<kappaFunction> operator()(Istream& is) const
        {
            //const word name(is);
            const dictionary parentDict(is);

            return autoPtr<kappaFunction>
            (
                kappaFunction::New(name_, parentDict, mesh_)
            );
        }
    };

    //- Return a reference to the selected fvOption model
    static autoPtr<kappaFunction> New
    (
        const word& name,
        const dictionary& parentDict,
        const fvMesh& mesh
    );

    //-
    virtual ~kappaFunction();

    //-
    virtual void update() = 0;

    //-
    virtual tmp<surfaceScalarField> kappa() = 0;

    //-
    const word& name() const;

    // I/O

    //-
    virtual void writeData (Ostream& ) const;

    //-
    virtual bool read (const dictionary& dict);
}; //class functionObject

}; //fv namespace

}; //Foam namespace

#endif

//END-OF-FILE
