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

#include "correctCentralACMIInterpolation.H"
#include "fvMesh.H"
#include "surfaceFields.H"
#include "cyclicACMIFvPatchFields.H"

void Foam::correctCentralACMIInterpolation(surfaceScalarField& nei_field)
{
    const fvMesh& mesh = nei_field.mesh();

    forAll(mesh.boundary(), iPatch)
    {
        if (isA<cyclicACMIFvPatch>(mesh.boundary()[iPatch]))
        {
            const cyclicACMIFvPatch& acmiPatch =
                refCast<const cyclicACMIFvPatch>(mesh.boundary()[iPatch]);

            if (acmiPatch.owner())
            {
                const cyclicACMIFvPatch& acmiOwnPatch = acmiPatch;
                const cyclicACMIFvPatch& acmiNeiPatch = acmiPatch.neighbPatch();

                label nonOverlapOwnID = acmiOwnPatch.nonOverlapPatchID();
                label nonOverlapNeiID = acmiNeiPatch.nonOverlapPatchID();

                //correct interpolation for nei field on own ACMI Patch
                nei_field.boundaryFieldRef()[acmiOwnPatch.index()] +=
                    (1.0 - acmiOwnPatch.AMI().srcWeightsSum())*nei_field.boundaryField()[nonOverlapOwnID];

                    //correct interpolation for nei field on nei ACMI Patch
                    //weights on nei ACMI Patch = acmiOwnPatch.AMI().tgtWeightsSum()
                    nei_field.boundaryFieldRef()[acmiNeiPatch.index()] +=
                        (1.0 - acmiOwnPatch.AMI().tgtWeightsSum())*nei_field.boundaryField()[nonOverlapNeiID];
            }
        }
    }
}

void Foam::correctCentralACMIInterpolation(surfaceVectorField& nei_field)
{
    const fvMesh& mesh = nei_field.mesh();

    forAll(mesh.boundary(), iPatch)
    {
        if (isA<cyclicACMIFvPatch>(mesh.boundary()[iPatch]))
        {
            const cyclicACMIFvPatch& acmiPatch =
                refCast<const cyclicACMIFvPatch>(mesh.boundary()[iPatch]);

            if (acmiPatch.owner())
            {
                const cyclicACMIFvPatch& acmiOwnPatch = acmiPatch;
                const cyclicACMIFvPatch& acmiNeiPatch = acmiPatch.neighbPatch();

                label nonOverlapOwnID = acmiOwnPatch.nonOverlapPatchID();
                label nonOverlapNeiID = acmiNeiPatch.nonOverlapPatchID();
                
                //correct interpolation for nei field on own ACMI Patch
                nei_field.boundaryFieldRef()[acmiOwnPatch.index()] +=
                    (1.0 - acmiOwnPatch.AMI().srcWeightsSum())*nei_field.boundaryField()[nonOverlapOwnID];

                //correct interpolation for nei field on nei ACMI Patch
                //weights on nei ACMI Patch = acmiOwnPatch.AMI().tgtWeightsSum()
                nei_field.boundaryFieldRef()[acmiNeiPatch.index()] +=
                    (1.0 - acmiOwnPatch.AMI().tgtWeightsSum())*nei_field.boundaryField()[nonOverlapNeiID];
            }
        }
    }
}

//
//END-OF-FILE
//
