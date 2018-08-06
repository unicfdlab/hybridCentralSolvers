{
    volScalarField& rDeltaT = trDeltaT.ref();

    if (maxCoPtr.empty())
    {
        maxCoPtr = Function1<scalar>::New("maxCo", runTime.controlDict());
    }
    maxCo = maxCoPtr->value(runTime.value());
    
    
    if (rDeltaTSmoothingCoeffPtr.empty())
    {
        rDeltaTSmoothingCoeffPtr = Function1<scalar>::New
        (
            "rDeltaTSmoothingCoeff",
            runTime.controlDict()
        );
    }
    scalar rDeltaTSmoothingCoeff
    (
        rDeltaTSmoothingCoeffPtr->value(runTime.value()) // 0.02 is default
    );
    
    
    if (rDeltaTDampingCoeffPtr.empty())
    {
        rDeltaTDampingCoeffPtr = Function1<scalar>::New
        (
            "rDeltaTDampingCoeff",
            runTime.controlDict()
        );
    }
    scalar rDeltaTDampingCoeff
    (
        rDeltaTDampingCoeffPtr->value(runTime.value()) // 1.0 is default - no damping
    );
    
    
    scalar maxDeltaT
    (
        readScalar(runTime.controlDict().lookup("maxDeltaT")) // GREAT - is a default
    );
    
    
    volScalarField rDeltaT0("rDeltaT0", rDeltaT);
    if (initialDeltaT < 0.0)
    {
        runTime.controlDict().lookup("initialDeltaT") >> initialDeltaT;
        rDeltaT.primitiveFieldRef() = 1.0 / initialDeltaT;
        rDeltaT0.primitiveFieldRef() = rDeltaT.primitiveField();
    }
    else
    {
        if ( CourantType == "cellCourant" )
        {
            // Set the reciprocal time-step from the local Courant number
    
            rDeltaT.ref() = max
            (
                1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
                (
                    0.5 * fvc::surfaceSum(mag(phi))().internalField() 
                    / mesh.V() / rho.internalField()
                ) / maxCo
            );
        }
        else if ( CourantType == "faceCourant" )
        {
            if (mesh.nInternalFaces())
            {
                surfaceScalarField phiOverRho =
                mag
                (
                    phi / (rho_own * alpha_own + rho_nei * alpha_nei)
                );
            
                surfaceScalarField SfUfbyDelta =
                    mesh.surfaceInterpolation::deltaCoeffs()*phiOverRho*2.0;
            
                forAll(mesh.V(), iCell)
                {
                    const labelList& cellFaces = mesh.cells()[iCell];
                    
                    scalar maxCellRDeltaT = -1.0;
                    scalar faceRDeltaT = maxCellRDeltaT;
                    label faceId = -1;
                    forAll(cellFaces, jFace)
                    {
                        faceId = cellFaces[jFace];
                        if (mesh.isInternalFace(faceId))
                        {
                            faceRDeltaT = SfUfbyDelta.internalField()[faceId] / mesh.magSf().internalField()[faceId] / maxCo;
                        }
                        if (faceRDeltaT > maxCellRDeltaT)
                        {
                            maxCellRDeltaT = faceRDeltaT;
                        }
                    }
                    rDeltaT.primitiveFieldRef()[iCell] = maxCellRDeltaT;
                }
                forAll(SfUfbyDelta.boundaryField(), kPatch)
                {
                    if(polyPatch::constraintType(mesh.boundaryMesh().types()[kPatch]))
                    {
                    }
                    else
                    {
                        const fvsPatchField<scalar>& pSfUfbyDelta = SfUfbyDelta.boundaryField()[kPatch];
                        const labelList& patchCells = mesh.boundary()[kPatch].faceCells();
                        scalar faceRDeltaT = -1;
                        forAll(pSfUfbyDelta, qFace)
                        {
                            faceRDeltaT = pSfUfbyDelta[qFace] / mesh.magSf().boundaryField()[kPatch][qFace] / maxCo;
                            if (faceRDeltaT > rDeltaT.primitiveField()[patchCells[qFace]])
                            {
                                rDeltaT.primitiveFieldRef()[patchCells[qFace]] = faceRDeltaT;
                            }
                        }
                    }
                }
           }
        
        }
        else
        {
            FatalErrorIn
            (
                "setRDeltaT.H:"
            )
            << "Wrong type of Courant criterion: " << CourantType
            << endl << " must be one of:" 
            << endl << "1) cellCourant"
            << endl << "2) faceCourant"
            << endl << abort(FatalError);
        }
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/(rDeltaT.primitiveField() + 1.0e-12))
        << ", " << gMax(1/(rDeltaT.primitiveField() + 1.0e-12)) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
    // && runTime.timeIndex() > runTime.startTimeIndex() + 1
     && runTime.timeIndex() > runTime.startTimeIndex() 
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}
