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
    
    Info<< "Constructing compressibleTwoPhaseMixtureThermo\n" << endl;
    compressibleTwoPhaseMixtureThermo thermo(mesh);
    
    volScalarField& YLiq = thermo.YLiq();
    volScalarField& YGas = thermo.YGas();
    volScalarField& YbarLiq = thermo.YbarLiq();
    volScalarField& YbarGas = thermo.YbarGas();
    
    Info<< "Reading thermophysical properties\n" << endl;
    volScalarField& p = thermo.p();
    volScalarField& T = thermo.T();
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

    volScalarField c
    (
        "cSound",
        sqrt(thermo.Cp()/thermo.Cv() / psi)
    );
    
    /*
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
    
    volScalarField EkChange
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
    */
    
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
    
    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
    
    DynamicList<label> badQualityCells(0);
    
    Info<< "Declare and initialise the cumulative continuity error\n" << endl;
    scalar cumulativeContErr = 0;
    
    const dimensionedScalar& pMin = thermo.pMin();
    
    Switch updateEnergyInPISO (false);
    word fluxScheme (word::null);
    Switch noFluxWeighting(false);

    volScalarField KChange
    (
        "KChange",
        rho*K / runTime.deltaT() * 0.0
    );
        
    volScalarField KSource
    (
        "KSource",
        p / runTime.deltaT() * 0.0
    );


