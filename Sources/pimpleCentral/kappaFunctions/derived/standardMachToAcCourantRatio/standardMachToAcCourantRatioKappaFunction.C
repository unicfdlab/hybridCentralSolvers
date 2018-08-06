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

#include "standardMachToAcCourantRatioKappaFunction.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvsPatchFields.H"
#include "correctCentralACMIInterpolation.H"
#include "localEulerDdtScheme.H"

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
    const surfaceScalarField& rho_own = mesh_.thisDb().lookupObject<surfaceScalarField>("rho_own");
    const surfaceScalarField& rho_nei = mesh_.thisDb().lookupObject<surfaceScalarField>("rho_nei");
    const surfaceScalarField& a_own   = mesh_.thisDb().lookupObject<surfaceScalarField>("alpha_own");
    const surfaceScalarField& a_nei   = mesh_.thisDb().lookupObject<surfaceScalarField>("alpha_nei");
    const surfaceScalarField& cf_own  = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_own");
    const surfaceScalarField& cf_nei  = mesh_.thisDb().lookupObject<surfaceScalarField>("cf_nei");
    const surfaceScalarField& uMagSf  = mesh_.thisDb().lookupObject<surfaceScalarField>("uMagSf");

    surfaceScalarField amaxSfbyDelta
    (
        mesh_.surfaceInterpolation::deltaCoeffs()*amaxSf
    );

    autoPtr<surfaceScalarField> FaceAcCourantPtr;

    if (fv::localEulerDdt::enabled(mesh_))
    {
        dimensionedScalar cDeltaT = runTime_.deltaT();

        cDeltaT.value() = 1.0 / gMax
            (
                mesh_.thisDb().lookupObject<volScalarField>(fv::localEulerDdt::rDeltaTName).internalField()
            );
//        surfaceScalarField rDeltaTf = linearInterpolate
//        (
//            mesh_.thisDb().lookupObject<volScalarField>(fv::localEulerDdt::rDeltaTName)
//        );

        FaceAcCourantPtr.set
        (
            new surfaceScalarField
            (
                "FaceAcCourant",
                //amaxSfbyDelta/uMagSf / rDeltaTf
                amaxSfbyDelta/uMagSf * cDeltaT
            )
        );
    }
    else
    {
        FaceAcCourantPtr.set
        (
            new surfaceScalarField
            (
                "FaceAcCourant",
                amaxSfbyDelta/uMagSf * runTime_.deltaT()
            )
        );
    }

    const surfaceScalarField& FaceAcCourant = FaceAcCourantPtr();

    Info << "max/min FaceAcCourant: " << max(FaceAcCourant).value() << "/" << min(FaceAcCourant).value() << endl;

    surfaceScalarField Maf
    (
        "Maf",
        max
        (
            (mag(phi) / (rho_own*a_own + rho_nei*a_nei) / uMagSf)
            / (cf_own*a_own + cf_nei*a_nei),
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
    surfaceScalarField& kappa = tKappa.ref();

    resetCoupledBoundaries(kappa);

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
