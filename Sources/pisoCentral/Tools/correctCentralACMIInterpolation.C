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


