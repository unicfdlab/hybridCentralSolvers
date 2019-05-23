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

#include "simplePhasePsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
//#include "fixedJumpFvPatchFields.H"
//#include "fixedJumpAMIFvPatchFields.H"
//#include "EnthalpyJumpFvPatchScalarField.H"
//#include "EnthalpyJumpAMIFvPatchScalarField.H"

void Foam::simplePhasePsiThermo::heBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hbf = h.boundaryFieldRef();

    forAll(hbf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hbf[patchi]).gradient()
                = hbf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hbf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hbf[patchi]).refGrad()
                = hbf[patchi].fvPatchField::snGrad();
        }
    }
}


Foam::wordList Foam::simplePhasePsiThermo::heBoundaryBaseTypes()
{
    const volScalarField::Boundary& tbf =
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
    const volScalarField::Boundary& tbf =
    this->T_.boundaryField();

    wordList hbt = tbf.types();

    forAll(tbf, patchi)
    {
	if (isA<fixedValueFvPatchScalarField>(tbf[patchi]))
	{
	    hbt[patchi] = fixedEnergyFvPatchScalarField::typeName;
	}
	else if
	(
	    isA<zeroGradientFvPatchScalarField>(tbf[patchi])
	    || isA<fixedGradientFvPatchScalarField>(tbf[patchi])
	)
	{
	    hbt[patchi] = gradientEnergyFvPatchScalarField::typeName;
	}
	else if (isA<mixedFvPatchScalarField>(tbf[patchi]))
	{
	    hbt[patchi] = mixedEnergyFvPatchScalarField::typeName;
	}
//	else if (isA<fixedJumpFvPatchScalarField>(tbf[patchi]))
//	{
//	    hbt[patchi] = EnergyJumpFvPatchScalarField::typeName;
//	}
//	else if (isA<fixedJumpAMIFvPatchScalarField>(tbf[patchi]))
//	{
//	    hbt[patchi] = EnergyJumpAMIFvPatchScalarField::typeName;
//	}
	else if (tbf[patchi].type() == "energyRegionCoupledFvPatchScalarField")
	{
	    hbt[patchi] = "energyRegionCoupledFvPatchScalarField";
	}
    }

    return hbt;
}

Foam::simplePhasePsiThermo::simplePhasePsiThermo(const fvMesh& mesh, const dictionary& dict)
:
    name_(dict.name()),
    mesh_(mesh),
    p_(mesh.thisDb().lookupObject<volScalarField>(dict.get<word>("p"))),
    T_(mesh.thisDb().lookupObject<volScalarField>(dict.get<word>("T"))),
    he_
    (
	IOobject
	(
	    ("he." + name_),
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
	    ("psi." + name_),
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
	    ("alpha." + name_),
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
	    ("rho." + name_),
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
	    ("rhoZero." + name_),
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
	he_.boundaryFieldRef()[patchi] ==
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
	fvPatchScalarField& ppsi     = psi_.boundaryFieldRef()[patchi];
	fvPatchScalarField& palpha   = alpha_.boundaryFieldRef()[patchi];
	fvPatchScalarField& phe      = he_.boundaryFieldRef()[patchi];
	fvPatchScalarField& prho     = rho_.boundaryFieldRef()[patchi];
	fvPatchScalarField& prhoZero = rhoZero_.boundaryFieldRef()[patchi];

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
		("he." + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    he_.dimensions()
	)
    );

    volScalarField& he = the.ref();
    scalarField& heCells = he.ref();
    const scalarField& TCells = T.v();

    forAll(heCells, celli)
    {
	heCells[celli] =
	Cp_ * TCells[celli];
    }

    forAll(he.boundaryField(), patchi)
    {
	fvPatchScalarField& hep = he.boundaryFieldRef()[patchi];
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

    scalarField& he = the.ref();

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
    scalarField& he = the.ref();

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

    volScalarField& cp = tCp.ref();

    forAll(this->T_, celli)
    {
	cp[celli] = Cp_;
    }

    forAll(this->T_.boundaryField(), patchi)
    {
	//const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pCp = cp.boundaryFieldRef()[patchi];

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
    scalarField& Cp = tCp.ref();

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
		("Cv." + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimEnergy/dimMass/dimTemperature
	)
    );

    volScalarField& cv = tCv.ref();

    forAll(this->T_, celli)
    {
	cv[celli] = Cv_;
    }

    forAll(this->T_.boundaryField(), patchi)
    {
	//const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pCv = cv.boundaryFieldRef()[patchi];

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
    scalarField& Cv = tCv.ref();

    forAll(T, facei)
    {
	Cv[facei] =
	    Cv_;
    }

    return tCv;
}

Foam::tmp<Foam::volScalarField> Foam::simplePhasePsiThermo::gamma() const
{

    const scalar locGamma = Cp_ / Cv_;

    tmp<volScalarField> tGamma
    (
	new volScalarField
	(
	    IOobject
	    (
		("gamma." + name_),
		mesh_.time().timeName(),
		mesh_,
		IOobject::NO_READ,
		IOobject::NO_WRITE
	    ),
	    mesh_,
	    dimless
	)
    );

    volScalarField& gamma = tGamma.ref();

    forAll(this->T_, celli)
    {
	gamma[celli] = locGamma;
    }

    forAll(this->T_.boundaryField(), patchi)
    {
	//const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
	const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
	fvPatchScalarField& pGamma = gamma.boundaryFieldRef()[patchi];

	forAll(pT, facei)
	{
	    pGamma[facei] = locGamma;
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

    const scalar locGamma = Cp_ / Cv_;

    tmp<scalarField> tGamma(new scalarField(T.size()));
    scalarField& gamma = tGamma.ref();

    forAll(T, facei)
    {
	gamma[facei] =
	    locGamma;
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
		("CpByCpv." + name_),
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
    tKappa.ref().rename("kappa." + name_);
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
    tKappaEff.ref().rename(name_ + "kappaEff");
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
    tAlphaEff.ref().rename("alphaEff." + name_);
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
