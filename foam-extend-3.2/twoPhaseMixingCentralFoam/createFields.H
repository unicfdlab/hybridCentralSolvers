    Info<< "Reading g\n" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));

    Info<< "Constructing compressibleTwoPhaseMixtureThermo\n" << endl;
    compressibleTwoPhaseMixtureThermo thermo(mesh);
    
    volScalarField& YLiq = thermo.YLiq();
    volScalarField& YGas = thermo.YGas();
    volScalarField& YbarLiq = thermo.YbarLiq();
    volScalarField& YbarGas = thermo.YbarGas();
    
    Info<< "Reading thermophysical properties\n" << endl;
    volScalarField& p = thermo.p();
    const volScalarField& T = thermo.T();
    volScalarField& psi = const_cast<volScalarField&>(thermo.psi());
    thermo.correct();
    
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    
    Info<< "Calculating mass flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho * U) & mesh.Sf()
    );
    
    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );
    
    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );
    
    volScalarField KChange
    (
        IOobject
        (
            "EkChange",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );
    
    Info<< "Creating field gradp\n" << endl;
    volVectorField gradp
    (
        IOobject
        (
            "gradp",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("gradp", p.dimensions()/dimLength, 0) * pTraits<vector>::zero
    );
    
    volVectorField limGradP
    (
        "limGradP",
        gradp
    );
    
    DynamicList<label> badQualityCells(0);
    
    volScalarField rhoHat
    (
        "rhoHat",
        rho - psi*p
    );
    
    const dimensionedScalar& pMin = thermo.pMin();
    
    scalar Cid
    (
        1.0 
        - 
        min(max(readScalar(mesh.solutionDict().subDict("PIMPLE").lookup("Cid")),0.0),1.0)
    );
    


