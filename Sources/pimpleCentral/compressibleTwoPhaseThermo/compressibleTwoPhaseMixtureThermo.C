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

#include "compressibleTwoPhaseMixtureThermo.H"
#include "fixedEnergyFvPatchScalarField.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleTwoPhaseMixtureThermo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::compressibleTwoPhaseMixtureThermo::heBoundaryCorrection(volScalarField& h)
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixtureThermo::compressibleTwoPhaseMixtureThermo
(
    const fvMesh& mesh
)
:
    rhoThermo(mesh, word::null),
    compressibleTwoPhaseMixture(mesh, *this),
    pMin_(this->lookup("pMin")),
    thermoLiq_(nullptr),
    thermoGas_(nullptr),
    he_
    (
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    //limit pressure before proceeding
    p_ = max(p_, pMin_);

    thermoLiq_ .reset
    (
        new simplePhasePsiThermo(mesh, this->subDict(liqPhaseName()))
    );

    thermoGas_.reset
    (
        new simplePhasePsiThermo(mesh, this->subDict(gasPhaseName()))
    );

    he_ = YLiq()*thermoLiq_->he() + YGas()*thermoGas_->he();

    heBoundaryCorrection(he_);

    correct();

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::compressibleTwoPhaseMixtureThermo::~compressibleTwoPhaseMixtureThermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibleTwoPhaseMixtureThermo::correct()
{
    //update temperature and heat capacities of liquid and gas

    const scalar epsilonH = 1.0e-5;

    forAll(T_, celli)
    {
        scalar CpLiq = thermoLiq_->Cp(p_[celli], T_[celli]);
        scalar CpGas = thermoGas_->Cp(p_[celli], T_[celli]);
        scalar Cpm =  CpLiq * YLiq()[celli] + CpGas * YGas()[celli];

        scalar deltaH = Cpm * T_[celli] - he_[celli];
        scalar F = deltaH * deltaH;
        scalar dFdT = 2.0 * Cpm * deltaH;

        while (mag(deltaH) >= epsilonH)
        {
            T_[celli] = T_[celli] - F / dFdT;

            CpLiq = thermoLiq_->Cp(p_[celli], T_[celli]);
            CpGas = thermoGas_->Cp(p_[celli], T_[celli]);
            Cpm =  CpLiq * YLiq()[celli] + CpGas * YGas()[celli];

            deltaH = Cpm * T_[celli] - he_[celli];
            F = deltaH * deltaH;
            dFdT = 2.0 * Cpm * deltaH;
        }
    }

    forAll(T_.boundaryField(), patchi)
    {
	fvPatchScalarField&       pT  = T_.boundaryFieldRef()[patchi];
	const fvPatchScalarField& pp  = p_.boundaryField()[patchi];
	fvPatchScalarField&      hep  = he_.boundaryFieldRef()[patchi];

	if (pT.fixesValue())
	{
	    hep == YLiq().boundaryField()[patchi] * thermoLiq_->he(pp, pT, patchi) +
		    YGas().boundaryField()[patchi] * thermoGas_->he(pp, pT, patchi);
	}
	else
	{
	    forAll(pT, facei)
	    {
		scalar CpLiq = thermoLiq_->Cp(pp[facei], pT[facei]);
		scalar CpGas = thermoGas_->Cp(pp[facei], pT[facei]);
		scalar Cpm =  CpLiq * YLiq().boundaryField()[patchi][facei] + CpGas * YGas().boundaryField()[patchi][facei];

                scalar deltaH = Cpm * pT[facei] - hep[facei];
		scalar F = deltaH * deltaH;
		scalar dFdT = 2.0 * Cpm * deltaH;
		while (mag(deltaH) >= epsilonH)
		{
		    pT[facei] = pT[facei] - F / dFdT;

		    CpLiq = thermoLiq_->Cp(pp[facei], pT[facei]);
		    CpGas = thermoGas_->Cp(pp[facei], pT[facei]);
		    Cpm =  CpLiq * YLiq().boundaryField()[patchi][facei] + CpGas * YGas().boundaryField()[patchi][facei];

		    deltaH = Cpm * pT[facei] - hep[facei];
		    F = deltaH * deltaH;
		    dFdT = 2.0 * Cpm * deltaH;
		}
	    }
	}
    }

    Info << "max/min T: " << max(T_).value() << "/" << min(T_).value() << endl;

    //correct properties of liquid
    thermoLiq_->correct();

    //correct properties of gas
    thermoGas_->correct();

    //correct mixture properties
    rhoEff() = 1.0 / (YLiq()/thermoLiq_->rho() + YGas()/thermoGas_->rho());

    updateVolFrac(thermoLiq_->rho(), thermoGas_->rho());

    const volScalarField& rhoGas = thermoGas_->rho();
    const volScalarField& rhoLiq = thermoLiq_->rho();
    const volScalarField& psiGas = thermoGas_->psi();
    const volScalarField& psiLiq = thermoLiq_->psi();

    volScalarField YLiqByRhoLiq = YLiq() / rhoLiq;
    volScalarField YGasByRhoGas = YGas() / rhoGas;

    psi_ =
        -Foam::pow(YLiqByRhoLiq + YGasByRhoGas,-2.0)*
        (
            - (YLiqByRhoLiq / rhoLiq) * psiLiq
            - (YGasByRhoGas / rhoGas) * psiGas
        );

    mu_ = YbarLiq()*thermoLiq_->mu() + YbarGas()*thermoGas_->mu();
    alpha_ = YbarLiq()*thermoLiq_->alpha() + YbarGas()*thermoGas_->alpha();

    //mu_ = YLiq()*thermoLiq_->mu() + YGas()*thermoGas_->mu();
    //alpha_ = YLiq()*thermoLiq_->alpha() + YGas()*thermoGas_->alpha();
}


Foam::word Foam::compressibleTwoPhaseMixtureThermo::thermoName() const
{
    return "compressibleTwoPhaseMixtureThermo";
}

bool Foam::compressibleTwoPhaseMixtureThermo::incompressible() const
{
    return false;
}


bool Foam::compressibleTwoPhaseMixtureThermo::isochoric() const
{
    return false;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return YLiq()*thermoLiq_->he(p, T) + YGas()*thermoGas_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    return
        scalarField(YLiq(), cells)*thermoLiq_->he(p, T, cells)
      + scalarField(YGas(), cells)*thermoGas_->he(p, T, cells);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->he(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->he(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::hc() const
{
    return YLiq()*thermoLiq_->hc() + YGas()*thermoGas_->hc();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    notImplemented("compressibleTwoPhaseMixtureThermo::THE(...)");
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    notImplemented("compressibleTwoPhaseMixtureThermo::THE(...)");
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::Cp() const
{
    return YLiq()*thermoLiq_->Cp() + YGas()*thermoGas_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->Cp(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->Cp(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::Cv() const
{
    return YLiq()*thermoLiq_->Cv() + YGas()*thermoGas_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->Cv(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->Cv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::gamma() const
{
    return YLiq()*thermoLiq_->gamma() + YGas()*thermoGas_->gamma();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->gamma(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->gamma(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::Cpv() const
{
    return YLiq()*thermoLiq_->Cpv() + YGas()*thermoGas_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->Cpv(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->Cpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::CpByCpv() const
{
    return
        YLiq()*thermoLiq_->CpByCpv()
      + YGas()*thermoGas_->CpByCpv();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->CpByCpv(p, T, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->CpByCpv(p, T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::kappa() const
{
    return YLiq()*thermoLiq_->kappa() + YGas()*thermoGas_->kappa();
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::kappa
(
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->kappa(patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return
        YLiq()*thermoLiq_->kappaEff(alphat)
      + YGas()*thermoGas_->kappaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->kappaEff(alphat, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->kappaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return
        YLiq()*thermoLiq_->alphaEff(alphat)
      + YGas()*thermoGas_->alphaEff(alphat);
}


Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        YLiq().boundaryField()[patchi]*thermoLiq_->alphaEff(alphat, patchi)
      + YGas().boundaryField()[patchi]*thermoGas_->alphaEff(alphat, patchi)
    ;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::mu() const
{
    return tmp<volScalarField>(this->mu_);
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField> Foam::compressibleTwoPhaseMixtureThermo::rho() const
{
    return rhoEff_;
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::rho(const label patchi) const
{
    return rhoEff_.boundaryField()[patchi];
}

const Foam::dimensionedScalar& Foam::compressibleTwoPhaseMixtureThermo::pMin() const
{
    return pMin_;
}

Foam::volScalarField& Foam::compressibleTwoPhaseMixtureThermo::h()
{
    return he_;
}

const Foam::volScalarField& Foam::compressibleTwoPhaseMixtureThermo::h() const
{
    return he_;
}

Foam::tmp<Foam::scalarField> Foam::compressibleTwoPhaseMixtureThermo::h
(
    const scalarField& T,
    const label patchi
) const
{
    return tmp<scalarField>(he_.boundaryField()[patchi]);
}

void Foam::compressibleTwoPhaseMixtureThermo::correctRealDensities()
{
    volScalarField pLimited
    (
        "pLim",
        max(pMin_, p_)
    );
    thermoLiq_->correctDensity(pLimited);
    thermoGas_->correctDensity(pLimited);
}


// ************************************************************************* //
