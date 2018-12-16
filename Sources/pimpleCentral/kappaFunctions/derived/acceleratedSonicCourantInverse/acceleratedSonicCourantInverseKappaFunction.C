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

#include "acceleratedSonicCourantInverseKappaFunction.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "localEulerDdtScheme.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(acceleratedSonicCourantInverseKappaFunction, 0);
    addToRunTimeSelectionTable
    (
        kappaFunction,
        acceleratedSonicCourantInverseKappaFunction,
        dictionary
    );
}
}

namespace Foam
{
namespace fv
{

acceleratedSonicCourantInverseKappaFunction::acceleratedSonicCourantInverseKappaFunction
(
    const word& name,
    const word& type,
    const dictionary& parentDict,
    const fvMesh& mesh
)
:
    kappaFunction (name, type, parentDict, mesh),
    power_(1.0),
    reinterpolateToFaces_(false)
{
    this->read(this->coeffs_);
}

acceleratedSonicCourantInverseKappaFunction::~acceleratedSonicCourantInverseKappaFunction()
{
}

void acceleratedSonicCourantInverseKappaFunction::update()
{
}

tmp<surfaceScalarField> acceleratedSonicCourantInverseKappaFunction::kappa()
{
    const surfaceScalarField& cf_own = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_own");
    const surfaceScalarField& cf_nei = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_nei");

    surfaceScalarField cfbyDelta
    (
        mesh_.surfaceInterpolation::deltaCoeffs()
        *
        (
            max(cf_own,cf_nei)
        )
    );

    dimensionedScalar cDeltaT = runTime_.deltaT();

    if (fv::localEulerDdt::enabled(mesh_))
    {
        cDeltaT.value() = 1.0 / gMax
        (
            mesh_.thisDb().lookupObject<volScalarField>(fv::localEulerDdt::rDeltaTName)
        );
    }

    surfaceScalarField FaceSonicCourant
    (
        "FaceSonicCourant",
        (cfbyDelta * cDeltaT)
    );

    Info << "max/min FaceSonicCourant: " << max(FaceSonicCourant).value() << "/" << min(FaceSonicCourant).value() << endl;

    tmp<surfaceScalarField> tKappa
    (
        pow
        (
            min(1.0 / FaceSonicCourant, 1.0),
            power_
        )
    );
    surfaceScalarField& kappa = tKappa.ref();
    kappa.setOriented(true);

    reinterpolateToFaces(kappa);

    resetCoupledBoundaries(kappa);

    writeMaxMinKappa(kappa);

    return tKappa;
}

void acceleratedSonicCourantInverseKappaFunction::reinterpolateToFaces(surfaceScalarField& kappa)
{
    if (reinterpolateToFaces_)
    {
        volScalarField minKappa
        (
            IOobject
            (
                "minKappa",
                runTime_.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("minkappa", dimless,  0.0)
        );

        forAll(mesh_.cells(), iCell)
        {
            scalar minKappaCell =  1.01;

            const labelList& cellFaces = mesh_.cells()[iCell];

            forAll(cellFaces, iFace)
            {
                if (mesh_.isInternalFace(cellFaces[iFace]))
                {
                    if (kappa[cellFaces[iFace]] < minKappaCell)
                    {
                        minKappaCell = kappa[cellFaces[iFace]];
                    }
                }
            }

            minKappa[iCell] = minKappaCell;
        }

        forAll(kappa.boundaryField(), iPatch)
        {
            forAll(kappa.boundaryField()[iPatch], iFace)
            {
                minKappa.boundaryFieldRef()[iPatch][iFace] =
                    kappa.boundaryField()[iPatch][iFace];
            }
        }

        kappa = linearInterpolate(minKappa);
    }
}

void acceleratedSonicCourantInverseKappaFunction::writeData (Ostream& os) const
{
    kappaFunction::writeData(os);
}

bool acceleratedSonicCourantInverseKappaFunction::read(const dictionary& dict)
{
    if (kappaFunction::read(dict))
    {
        dict.lookup("power") >> power_;

        reinterpolateToFaces_ = dict.lookupOrDefault<Switch>("reinterpolateToFaces", false);

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
