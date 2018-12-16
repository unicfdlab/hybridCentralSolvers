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

#include "sonicCourantInverseKappaFunction.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvsPatchFields.H"
#include "localEulerDdtScheme.H"

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
    const surfaceScalarField& cf_own = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_own");
    const surfaceScalarField& cf_nei = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_nei");
    const surfaceScalarField& a_own   = mesh_.thisDb().lookupObject<surfaceScalarField>("alpha_own");
    const surfaceScalarField& a_nei   = mesh_.thisDb().lookupObject<surfaceScalarField>("alpha_nei");


    surfaceScalarField cfbyDelta
    (
        mesh_.surfaceInterpolation::deltaCoeffs()
        *
        (
            cf_own*a_own + cf_nei*a_nei
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
        min
        (
            1.0 / FaceSonicCourant,
            scalar(1.0)
        )
    );
    surfaceScalarField& kappa = tKappa.ref();
    kappa.setOriented(true);

    resetCoupledBoundaries(kappa);

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
