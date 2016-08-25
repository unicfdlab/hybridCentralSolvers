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

#include "standardMachToAcCourantRatioKappaFunction.H"
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
    defineTypeNameAndDebug(standardMachToAcCourantRatioKappaFunction, 0);
    addToRunTimeSelectionTable
    (
        kappaFunction,
        standardMachToAcCourantRatioKappaFunction,
        dictionary
    );
}
}

namespace Foam
{
namespace fv
{

standardMachToAcCourantRatioKappaFunction::standardMachToAcCourantRatioKappaFunction
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

standardMachToAcCourantRatioKappaFunction::~standardMachToAcCourantRatioKappaFunction()
{
}

void standardMachToAcCourantRatioKappaFunction::update()
{
}

tmp<surfaceScalarField> standardMachToAcCourantRatioKappaFunction::kappa()
{
    const surfaceScalarField& amaxSf  = mesh_.thisDb().lookupObject<surfaceScalarField>("amaxSf");
    const surfaceScalarField& phi     = mesh_.thisDb().lookupObject<surfaceScalarField>("phi");
    const surfaceScalarField& rho_pos = mesh_.thisDb().lookupObject<surfaceScalarField>("rho_pos");
    const surfaceScalarField& rho_neg = mesh_.thisDb().lookupObject<surfaceScalarField>("rho_neg");
    const surfaceScalarField& a_pos   = mesh_.thisDb().lookupObject<surfaceScalarField>("a_pos");
    const surfaceScalarField& a_neg   = mesh_.thisDb().lookupObject<surfaceScalarField>("a_neg");
    const surfaceScalarField& cSf_pos = mesh_.thisDb().lookupObject<surfaceScalarField>("cSf_pos");
    const surfaceScalarField& cSf_neg = mesh_.thisDb().lookupObject<surfaceScalarField>("cSf_neg");
    
    surfaceScalarField amaxSfbyDelta
    (
        mesh_.surfaceInterpolation::deltaCoeffs()*amaxSf
    );
    
    surfaceScalarField FaceAcCourant
    (
        "FaceAcCourant",
        (amaxSfbyDelta/mesh_.magSf() * runTime_.deltaT())
    );
    
    Info << "max/min FaceAcCourant: " << max(FaceAcCourant).value() << "/" << min(FaceAcCourant).value() << endl;
    
    surfaceScalarField Maf
    (
        max
        (
            mag(phi) / (rho_pos*a_pos + rho_neg*a_neg)
            / (cSf_pos*a_pos + cSf_neg*a_neg),
            scalar(0)
        )
    );
    
    Info << "max/min Maf: " << max(Maf).value() << "/" << min(Maf).value() << endl;
    
    tmp<surfaceScalarField> tKappa
    (
        min
        (
            Maf / FaceAcCourant,
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

void standardMachToAcCourantRatioKappaFunction::writeData (Ostream& os) const
{
    kappaFunction::writeData(os);
}

bool standardMachToAcCourantRatioKappaFunction::read(const dictionary& dict)
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

