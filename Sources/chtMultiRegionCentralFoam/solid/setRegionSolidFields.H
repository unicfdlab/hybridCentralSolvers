fvMesh& solidMesh = solidRegions[i];
solidThermo& solidThermo = thermos[i];

tmp<volScalarField> tSolidRho = solidThermo.rho();
const volScalarField& solidRho = tSolidRho();

tmp<volScalarField> tSolidCp = solidThermo.Cp();
const volScalarField& solidCp = tSolidCp();

tmp<volSymmTensorField> taniAlpha;
if (!solidThermo.isotropic())
{
    volSymmTensorField& aniAlpha = aniAlphas[i];
    tmp<volVectorField> tkappaByCp = solidThermo.Kappa()/solidCp;
    const coordinateSystem& coodSys = coordinates[i];

    aniAlpha.primitiveFieldRef() =
        coodSys.transformPrincipal
        (
           solidRegions[i].cellCentres(),
           tkappaByCp()
        );
    aniAlpha.correctBoundaryConditions();

    taniAlpha = tmp<volSymmTensorField>
    (
        new volSymmTensorField(aniAlpha)
    );
}


volScalarField& solidH = solidThermo.he();

const volScalarField& solidBetav = betavSolid[i];

fv::options& solidFvOptions = solidHeatSources[i];

