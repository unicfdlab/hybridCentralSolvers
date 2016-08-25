/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "kappaFunction.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "foamTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(kappaFunction, 0);
    defineRunTimeSelectionTable(kappaFunction, dictionary);
}
}

namespace Foam
{
namespace fv
{

autoPtr<kappaFunction> kappaFunction::New
(
    const word& name,
    const dictionary& parentDict,
    const fvMesh& mesh
)
{

    word modelType(parentDict.lookup("type"));
    
    Info<< "Selecting finite volume kappaFunction type " << modelType << endl;
    
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "kappaFunction::New(const word&, const dictionary&, const fvMesh&)"
        )   << "Unknown Model type " << modelType << nl << nl
        << "Valid model types are:" << nl
        << dictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }

    return autoPtr<kappaFunction>(cstrIter()(name, modelType, parentDict, mesh));
}

kappaFunction::kappaFunction(const word& name, const word& type, const dictionary& parentDict, const fvMesh& mesh)
:
    name_(name),
    scType_(type),
    mesh_(mesh),
    runTime_(mesh.time()),
    coeffs_(parentDict.subDict(scType_ + "Coeffs")),
    writeMaxMin_(false)
{
}

kappaFunction::~kappaFunction()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kappaFunction::writeMaxMinKappa (const surfaceScalarField& kappa)
{
    if (writeMaxMin_ && runTime_.outputTime())
    {
        volScalarField maxKappa
        (
            IOobject
            (
                "maxKappa",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("maxkappa", dimless,  0.0)
        );
        
        volScalarField minKappa
        (
            "minKappa",
            maxKappa*0.0
        );
        
        forAll(mesh_.cells(), iCell)
        {
            scalar maxKappaCell =  -0.01;
            scalar minKappaCell =  1.01;
            
            const labelList& cellFaces = mesh_.cells()[iCell];
            
            forAll(cellFaces, iFace)
            {
                if (mesh_.isInternalFace(cellFaces[iFace]))
                {
                    if (kappa[cellFaces[iFace]] > maxKappaCell)
                    {
                        maxKappaCell = kappa[cellFaces[iFace]];
                    }
                    if (kappa[cellFaces[iFace]] < minKappaCell)
                    {
                        minKappaCell = kappa[cellFaces[iFace]];
                    }
                }
            }
            
            maxKappa[iCell] = maxKappaCell;
            minKappa[iCell] = minKappaCell;
        }
        
        maxKappa.write();
        minKappa.write();
    }
}

void kappaFunction::writeData (Ostream& os) const
{
}

bool kappaFunction::read(const dictionary& dict)
{
    
    coeffs_.lookup("writeMaxMin") >> writeMaxMin_;
    
    return true;
}


const word& kappaFunction::name() const
{
    return name_;
}

}; //namespace fv

}; //namespace Foam


//END-OF-FILE

