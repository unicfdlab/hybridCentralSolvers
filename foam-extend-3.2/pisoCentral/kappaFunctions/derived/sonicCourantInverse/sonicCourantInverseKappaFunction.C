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

#include "sonicCourantInverseKappaFunction.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvsPatchFields.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(sonicCourantInverseKappaFunction, 0);
    addToRunTimeSelectionTable
    (
        kappaFunction,
        sonicCourantInverseKappaFunction,
        dictionary
    );
}
}

namespace Foam
{
namespace fv
{

sonicCourantInverseKappaFunction::sonicCourantInverseKappaFunction
(
    const word& name,
    const word& type,
    const dictionary& parentDict,
    const fvMesh& mesh
)
:
    kappaFunction (name, type, parentDict, mesh)
{
    this->read(this->coeffs_);
}

sonicCourantInverseKappaFunction::~sonicCourantInverseKappaFunction()
{
}

void sonicCourantInverseKappaFunction::update()
{
}

tmp<surfaceScalarField> sonicCourantInverseKappaFunction::kappa()
{
    const surfaceScalarField& cSf_pos = mesh_.thisDb().lookupObject<surfaceScalarField>("cSf_pos");
    const surfaceScalarField& cSf_neg = mesh_.thisDb().lookupObject<surfaceScalarField>("cSf_neg");
    
    surfaceScalarField cSfbyDelta
    (
        mesh_.surfaceInterpolation::deltaCoeffs()
        *
        (
            max(cSf_pos,cSf_neg)
        )
    );
    
    surfaceScalarField FaceSonicCourant
    (
        "FaceSonicCourant",
        (cSfbyDelta/mesh_.magSf() * runTime_.deltaT())
    );
    
    Info << "max/min FaceSonicCourant: " << max(FaceSonicCourant).value() << "/" << min(FaceSonicCourant).value() << endl;
    
    tmp<surfaceScalarField> tKappa
    (
        min
        (
            1.0 / FaceSonicCourant,
            scalar(1.0)
        )
    );
    surfaceScalarField& kappa = tKappa();

    forAll(kappa.boundaryField(), iPatch)
    {
        fvsPatchField<scalar>& kappapf = kappa.boundaryField()[iPatch];
        if (isA<coupledFvsPatchField<scalar> > (kappapf))
        {
            forAll (kappapf, iFace)
            {
                kappapf[iFace] = 0.0;
            }
        }
    }
    
    writeMaxMinKappa(kappa);

    return tKappa;
}

void sonicCourantInverseKappaFunction::writeData (Ostream& os) const
{
    kappaFunction::writeData(os);
}

bool sonicCourantInverseKappaFunction::read(const dictionary& dict)
{
    if (kappaFunction::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
    
    return true;
}

}; //namespace fv

}; //namespace Foam


//END-OF-FILE

