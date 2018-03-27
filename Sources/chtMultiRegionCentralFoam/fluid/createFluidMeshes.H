    const wordList fluidNames(rp["fluid"]);

    autoPtr<fvMesh> fluidMeshPtr;
    
    if (fluidNames.size() < 1)
    {
        FatalErrorIn
        (
            "chtMultiRegionCentralFoam"
        )   << "No fluid regions in case "
            << "Add fluid regions" << nl
            << exit(FatalError);

    }
    else
    {
        Info<< "Create fluid mesh for region " << fluidNames[0]
            << " for time = " << runTime.timeName() << nl << endl;

        fluidMeshPtr.reset
        (
            new fvMesh
            (
                IOobject
                (
                    fluidNames[0],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
    
    fvMesh& fluidMesh = fluidMeshPtr();
    fvMesh& mesh = fluidMesh;


