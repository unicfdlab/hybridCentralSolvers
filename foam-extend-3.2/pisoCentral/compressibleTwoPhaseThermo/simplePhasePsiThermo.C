#include "simplePhasePsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnthalpyFvPatchScalarField.H"
#include "gradientEnthalpyFvPatchScalarField.H"
#include "mixedEnthalpyFvPatchScalarField.H"
//#include "fixedJumpFvPatchFields.H"
//#include "fixedJumpAMIFvPatchFields.H"
//#include "EnthalpyJumpFvPatchScalarField.H"
//#include "EnthalpyJumpAMIFvPatchScalarField.H"

void Foam::simplePhasePsiThermo::heBoundaryCorrection(volScalarField& h)
{
    volScalarField::GeometricBoundaryField& hbf = h.boundaryField();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnthalpyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnthalpyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnthalpyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


Foam::wordList Foam::simplePhasePsiThermo::heBoundaryBaseTypes()
{
    const volScalarField::GeometricBoundaryField& tbf =
    this->T_.boundaryField();

    wordList hbt(tbf.size(), word::null);

    forAll(tbf, patchi)
    {
//	if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
//	{
//	    const fixedJumpFvPatchScalarField& pf =
//		dynamic_cast<const fixedJumpFvPatchScalarField&>(tbf[patchi]);
//	    
//	    hbt[patchi] = pf.interfaceFieldType();
//	}
//	else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
//	{
//	    const fixedJumpAMIFvPatchScalarField& pf =
//	    dynamic_cast<const fixedJumpAMIFvPatchScalarField&>
//	    (
//		tbf[patchi]
//	    );
//	    
//	    hbt[patchi] = pf.interfaceFieldType();
//	}
    }
    
    return hbt;
}

Foam::wordList Foam::simplePhasePsiThermo::heBoundaryTypes()
{
    const volScalarField::GeometricBoundaryField& tbf =
    this->T_.boundaryField();
    
    wordList hbt = tbf.types();
    
    forAll(tbf, patchi)
    {
	if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
	{
	    hbt[patchi] = fixedEnthalpyFvPatchScalarField::typeName;
	}
	else if
	(
	    isA<zeroGradientFvPatchScalarField>(tbf[patchi])
	    || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
	)
	{
	    hbt[patchi] = gradientEnthalpyFvPatchScalarField::typeName;
	}
	else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
	{
	    hbt[patchi] = mixedEnthalpyFvPatchScalarField::typeName;
	}
//	else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
//	{
//	    hbt[patchi] = EnthalpyJumpFvPatchScalarField::typeName;
//	}
//	else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
//	{
//	    hbt[patchi] = EnthalpyJumpAMIFvPatchScalarField::typeName;
//	}
	else if (tbf[patchi].type() == "enthalpyRegionCoupledFvPatchScalarField")
	{
	    hbt[patchi] = "enthalpyRegionCoupledFvPatchScalarField";
	}
    }
    
    return hbt;
}

Foam::simplePhasePsiThermo::simplePhasePsiThermo(const fvMesh& mesh, const dictionary& dict)
:
    name_(dict.name()),
    mesh_(mesh),
    p_(mesh.thisDb().lookupObject<volScalarField>(dict.lookup("p"))),
    T_(mesh.thisDb().lookupObject<volScalarField>(dict.lookup("T"))),
    he_
    (
	IOobject
	(
	    ("he" + name_),
	    mesh.time().timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimEnergy/dimMass,
	this->heBoundaryTypes()
    ),
    psi_
    (
	IOobject
	(
	    ("psi" + name_),
	    mesh.time().timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimTime*dimTime / dimLength / dimLength
    ),
    alpha_
    (
	IOobject
	(
	    ("alpha" + name_),
	    mesh.time().timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimMass / dimTime / dimLength
    ),
    rho_
    (
	IOobject
	(
	    ("rho" + name_),
	    mesh.time().timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimMass / (dimLength*dimLength*dimLength)
    ),
    rhoZero_
    (
	IOobject
	(
	    ("rhoZero" + name_),
	    mesh.time().timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimMass / (dimLength*dimLength*dimLength)
    )
{
    R_   = readScalar(dict.lookup("R"));
    Cp_  = readScalar(dict.lookup("Cp"));
    Cv_  = readScalar(dict.lookup("Cv"));
    mu_  = readScalar(dict.lookup("mu"));
    Pr_  = readScalar(dict.lookup("Pr"));
    rho0_= readScalar(dict.lookup("rho0"));
    p0_  = readScalar(dict.lookup("p0"));
    
    this->he_.operator= (this->he(p_,T_)());
    
    forAll(he_.boundaryField(), patchi)
    {
	he_.boundaryField()[patchi] == 
	    this->he
	    (
		p_.boundaryField()[patchi],
		T_.boundaryField()[patchi],
		patchi
	    );
    }
    
    heBoundaryCorrection(he_);
    
    correct();
    
}

Foam::simplePhasePsiThermo::~simplePhasePsiThermo()
{
}

void Foam::simplePhasePsiThermo::correct()
{
    forAll(T_, celli)
    {
	psi_[celli] = 1.0 / (R_ * T_[celli]);
	alpha_[celli] = mu_ / Pr_;
	rhoZero_[celli] = rho0_ - psi_[celli] * p0_;
	rho_[celli] = rhoZero_[celli] + psi_[celli] * p_[celli];
    }
    
    forAll(T_.boundaryField(), patchi)
    {
	const fvPatchScalarField& pT = T_.boundaryField()[patchi];
	const fvPatchScalarField& pp = p_.boundaryField()[patchi];
	fvPatchScalarField& ppsi     = psi_.boundaryField()[patchi];
	fvPatchScalarField& palpha   = alpha_.boundaryField()[patchi];
	fvPatchScalarField& phe      = he_.boundaryField()[patchi];
	fvPatchScalarField& prho     = rho_.boundaryField()[patchi];
	fvPatchScalarField& prhoZero = rhoZero_.boundaryField()[patchi];
	
	if (pT.fixesValue())
	{
	    forAll(pT, facei)
	    {
		phe[facei] = pT[facei] * Cp_;
	    }
	}
	forAll(pT, facei)
	{
	    ppsi[facei] = 1.0 / (R_ * pT[facei]);
	    palpha[facei] = mu_ / Pr_;
	    prhoZero[facei] = rho0_ - ppsi[facei] * p0_;
	    prho[facei] = prhoZero[facei] + pp[facei] * ppsi[facei];
	}
    }
}


Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> the
    (
	new volScalarField
	(
	    IOobject
	    (
		("he" + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    he_.dimensions()
	)
    );

    volScalarField& he = the();
    scalarField& heCells = he.internalField();
    const scalarField& TCells = T.internalField();

    forAll(heCells, celli)
    {
	heCells[celli] =
	Cp_ * TCells[celli];
    }

    forAll(he.boundaryField(), patchi)
    {
	fvPatchScalarField& hep = he.boundaryField()[patchi];
	const scalarField&  Tp  = T_.boundaryField()[patchi];
	
	forAll(hep, facei)
	{
	    hep[facei] =
	     Cp_ * Tp[facei];
	}
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    
    scalarField& he = the();
    
    forAll(T, celli)
    {
	he[celli] = Cp_* T[celli];
    }
    
    return the;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the();

    forAll(T, facei)
    {
	he[facei] =
	    Cp_ * T[facei];
    }

    return the;
}

Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::hc() const
{
    notImplemented("simplePhasePsiThermo::hc(...)");
    return tmp<volScalarField> ( volScalarField::null() );
}

Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    notImplemented("simplePhasePsiThermo::THE(...)");
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    notImplemented("simplePhasePsiThermo::THE(...)");
    return T0;
}

Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::Cp() const
{
    tmp<volScalarField> tCp
    (
	new volScalarField
	(
	    IOobject
	    (
		("Cp" + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimEnergy/dimMass/dimTemperature
	)
    );
    
    volScalarField& cp = tCp();
    
    forAll(this->T_, celli)
    {
	cp[celli] = Cp_;
    }
    
    forAll(this->T_.boundaryField(), patchi)
    {
	//const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pCp = cp.boundaryField()[patchi];
	
	forAll(pT, facei)
	{
	    pCp[facei] = Cp_;
	}
    }
    
    return tCp;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCp(new scalarField(T.size()));
    scalarField& Cp = tCp();

    forAll(T, facei)
    {
	Cp[facei] =
	    Cp_;
    }

    return tCp;
}

Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::Cv() const
{
    tmp<volScalarField> tCv
    (
	new volScalarField
	(
	    IOobject
	    (
		("Cv" + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimEnergy/dimMass/dimTemperature
	)
    );
    
    volScalarField& cv = tCv();
    
    forAll(this->T_, celli)
    {
	cv[celli] = Cv_;
    }
    
    forAll(this->T_.boundaryField(), patchi)
    {
	//const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pCv = cv.boundaryField()[patchi];
	
	forAll(pT, facei)
	{
	    pCv[facei] = Cv_;
	}
    }
    
    return tCv;
}

Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCv(new scalarField(T.size()));
    scalarField& Cv = tCv();

    forAll(T, facei)
    {
	Cv[facei] =
	    Cv_;
    }

    return tCv;
}

Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::gamma() const
{
    tmp<volScalarField> tGamma
    (
	new volScalarField
	(
	    IOobject
	    (
		("gamma" + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimless
	)
    );
    
    volScalarField& gamma = tGamma();
    
    forAll(this->T_, celli)
    {
	gamma[celli] = Cp_ / Cv_;
    }
    
    forAll(this->T_.boundaryField(), patchi)
    {
	//const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pGamma = gamma.boundaryField()[patchi];
	
	forAll(pT, facei)
	{
	    pGamma[facei] = Cp_ / Cv_;
	}
    }
    
    return tGamma;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tGamma(new scalarField(T.size()));
    scalarField& gamma = tGamma();

    forAll(T, facei)
    {
	gamma[facei] =
	    Cp_ / Cv_ ;
    }

    return tGamma;
}

Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::Cpv() const
{
    return Cp();
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return Cp(p,T,patchi);
}


Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::CpByCpv() const
{
    tmp<volScalarField> tCpByCpv
    (
	new volScalarField
	(
	    IOobject
	    (
		("CpByCpv" + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimless
	)
    );
    
    notImplemented("simplePhasePsiThermo::CpByCpv()");
    
    return tCpByCpv;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tCpByCpv(new scalarField(T.size()));
    
    notImplemented("simplePhasePsiThermo::CpByCpv(...)");

    return tCpByCpv;
}

Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::kappa() const
{
    tmp<volScalarField> tKappa(Cp() * this->alpha_);
    tKappa().rename(name_ + "kappa");
    return tKappa;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::kappa
(
    const label patchi
) const
{
    return
	Cp
	(
	    this->p_.boundaryField()[patchi],
	    this->T_.boundaryField()[patchi],
	    patchi
	) * alpha_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<volScalarField> tKappaEff(Cp() * (this->alpha_ + alphat));
    tKappaEff().rename(name_ + "kappaEff");
    return tKappaEff;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
	Cp
	(
	    this->p_.boundaryField()[patchi],
	    this->T_.boundaryField()[patchi],
	    patchi
	) * (alpha_.boundaryField()[patchi] + alphat);
}


Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<volScalarField> tAlphaEff(this->alpha_ + alphat);
    tAlphaEff().rename(name_ + "alphaEff");
    return tAlphaEff;
}


Foam::tmp<Foam::scalarField> Foam::simplePhasePsiThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return (alpha_.boundaryField()[patchi] + alphat);
}


//
//END-OF-FILE
//

