# Contents

1. [Available solvers with hybrid approximation](#Available-solvers-with-hybrid-approximation)
2. [Meeting points for users and developers](#Meeting-points-for-users-and-developers)
3. [Available OpenFOAM versions](#Available-OpenFOAM-versions)
4. [Derived projects](#Derived-projects)
5. [Research studies where the library was useful](#Research-studies-where-the-library-was-useful)
6. [For citation](#For-citation)

# Available solvers with hybrid approximation
[To the contents](#Contents)

United collection of hybrid  Central solvers - one-phase, two-phase and multicomponent versions.

Only OpenFOAM+ version of the OpenFOAM technology is supported since 2018. Use branches digitef-dev-YYMM, where YYMM - corresponds to OpenFOAM+ version, for example 1812

1. **pimpleCentralFoam** - Pressure-based semi implicit compressible flow of perfect gas solver based on 
central-upwind schemes of Kurganov and Tadmor with LTS support for steady-state calculations

2. **rhoPimpleCentralFoam** - Pressure-based semi implicit compressible flow  of real gas solver based on 
central-upwind schemes of Kurganov and Tadmor and LTS support for steady-state calculations

3. **pimpleCentralDyMFoam** - Pressure-based semi implicit compressible flow of perfect gas solver based 
on central-upwind schemes of Kurganov and Tadmor with mesh motion and LTS support for steady-state calculations

4. **reactingPimpleCentralFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of 
Kurganov and Tadmor for combustion with chemical reactions and LTS support for steady-state calculations

5. **twoPhaseMixingCentralFoam** - Transient Eulerian two-phase solver. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface is not accounted.

6. **twoPhaseMixingCentralDyMFoam** - Transient Eulerian two-phase solver with dynamic meshes. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface is not accounted.

7. **chtMultiRegionCentralFoam** -     Pressure-based semi implicit solver, based on hybrid central-upwind schemes
    of Kurganov and Tadmor for conjugate simulation of compressible flows of perfect gas (Mach 
    number is ranging from 0 to 6) and solid body heat transfer.

8. **interTwoPhaseCentralFoam** - pressure-based solver for compressible 0-4 Mach numbers flow of two-phase media with account to viscosity and gravity. The solver utilizes VoF method for resolution of phase interface and ACID to calculate properties in the region where both phases are present. The solver is located in [exper/2phase](https://github.com/unicfdlab/hybridCentralSolvers/tree/exper/2phase) branch (OpenFOAM+ v. 2012).


# Meeting points for users and developers
[To the contents](#Contents)

You can discuss questions of [hybridCentralSolvers](https://github.com/unicfdlab/hybridCentralSolvers) usage at Telegram Group: https://t.me/hybridCentralSolvers 

There is a [ResearchGate](https://www.researchgate.net/) project dedicated to the development of [hybridCentralSolvers  library](https://www.researchgate.net/project/Development-and-implementation-of-hybrid-Density-Pressure-scheme-for-compressible-flows-simulation-in-OpenFOAM)

# Available OpenFOAM versions
[To the contents](#Contents)

The library is available for next version of OpenFOAM:
* OpenFOAM 3.1 - [master branch](https://github.com/unicfdlab/hybridCentralSolvers/tree/master)
* OpenFOAM 4.1 - [dev-of4.1 branch](https://github.com/unicfdlab/hybridCentralSolvers/tree/dev-of4.1)
* OpenFOAM 6   - [dev-of6 branch](https://github.com/unicfdlab/hybridCentralSolvers/tree/dev-of6)
* OpenFOAM+ 1812 - [digitef-dev-1812](https://github.com/unicfdlab/hybridCentralSolvers/tree/digitef-dev-1812)
* OpenFOAM+ 1912 - [digitef-dev-1912](https://github.com/unicfdlab/hybridCentralSolvers/tree/digitef-dev-1912)
* OpenFOAM+ 2012 - [digitef-dev-2012](https://github.com/unicfdlab/hybridCentralSolvers/tree/digitef-dev-2012)

# Derived projects
[To the contents](#Contents)

The library or approach were used in next projects:
* [multiRegionRectingPimpleCentralFoam](https://github.com/TonkomoLLC/hybridCentralSolvers/tree/master/OpenFOAM-4.1/multiRegionReactingPimpleCentralFoam) - the solver for coupled simulation of gas dynamics and heat transfer using hybrid KT/PIMPLE approximation of convective fluxes
* [adjointReactingRhoPimpleCentralFoam](https://github.com/clapointe2011/public/tree/master/discreteAdjointOpenFOAM/applications/solvers/adjoint/adjointReactingRhoPimpleCentralFoam) - the solver for adjoint shape optimization of region with gas flow modelled using hybrid KT/PIMPLE approximation of convective fluxes


# Research studies where the library was useful
[To the contents](#Contents)

If you want to see your research in this list, please write to [Issues](https://github.com/unicfdlab/hybridCentralSolvers/issues) .

## <p align="center"> >>>>> 2022 <<<<<  </p> 

| Title | Description |
|------|-------------|
|[Analysis of the ignition induced by shock wave focusing equipped with conical and hemispherical reflectors](https://www.researchgate.net/publication/355077511_Analysis_of_the_ignition_induced_by_shock_wave_focusing_equipped_with_conical_and_hemispherical_reflectors):  **Article** |![Temperature distribution](https://ars.els-cdn.com/content/image/1-s2.0-S001021802100506X-gr5.jpg)| 


## <p align="center"> >>>>> 2021 <<<<<  </p> 

| Title | Description |
|------|-------------|
|[The Eulerian–Lagrangian Approach for the Numerical Investigation of an Acoustic Field Generated by a High-Speed Gas-Droplet Flow](https://www.mdpi.com/2311-5521/6/8/274):  **Article**| ![Jet with particles Logo](https://www.mdpi.com/fluids/fluids-06-00274/article_deploy/html/images/fluids-06-00274-ag-550.jpg) |
|[Numerical Study of Thermoacoustic Waves in a Cavity under Rapid Wall Heating](https://www.researchgate.net/publication/355012423_Numerical_Study_of_Thermoacoustic_Waves_in_a_Cavity_under_Rapid_Wall_Heating):  **Article** |![Variation of the component of velocity ](https://media.springernature.com/lw685/springer-static/image/art%3A10.1134%2FS1995080221090122/MediaObjects/12202_2021_6494_Fig3_HTML.png?as=webp)|
|[Real-Gas Effects and Single-Phase Instabilities during Injection, Mixing and Combustion under High-Pressure Conditions](https://www.researchgate.net/publication/354224374_Real-Gas_Effects_and_Single-Phase_Instabilities_during_Injection_Mixing_and_Combustion_under_High-Pressure_Conditions):  **PhD Thesis** |![p-T diagram of a jet injection](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/p-T-jet.png)|
|[Dynamics of detonation transmission and propagation in a curved chamber: a numerical and experimental analysis](https://doi.org/10.1016/j.combustflame.2020.09.032):  **Article** |![Experiment vs calculation](https://ars.els-cdn.com/content/image/1-s2.0-S0010218020304168-gr2.jpg)|
|[Modelling of Supersonic and Subsonic Flows Using Hybrid PressureBased Solver in Openfoam](https://doi.org/10.11159/ffhmt21.107):  **Article** |![Schematic of Bluff body burner](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Bluff-body-sketch.png)|

## <p align="center"> >>>>> 2020 <<<<< </p>

| Title | Description |
|------|-------------|
|[Numerical simulation of transpiration cooling experiments in supersonic flow using OpenFOAM](https://link.springer.com/article/10.1007/s12567-019-00292-6) :  **Article** |![Schematic illustration of the applied porous interface model](https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs12567-019-00292-6/MediaObjects/12567_2019_292_Fig1_HTML.png?as=webp)|
|[On uncertainty quantification via the ensemble of independent numerical solutions](https://doi.org/10.1016/j.jocs.2020.101114):  **Article**  |![Flow scheme](https://ars.els-cdn.com/content/image/1-s2.0-S1877750319310695-gr1.jpg)|
|[Influence of hydrogen equivalence ratios on supersonic combustion based on large eddy simulations](https://doi.org/10.1016/j.ijhydene.2020.02.054):  **Article**  |![Model scramjet](https://ars.els-cdn.com/content/image/1-s2.0-S036031992030584X-gr1.jpg)|
|[A quasi-direct numerical simulation solver for compressible reacting flows](https://doi.org/10.1016/j.compfluid.2020.104718): **Article**|![Data exchange between OpenFOAM and Cantera](https://ars.els-cdn.com/content/image/1-s2.0-S0045793020302887-gr1.jpg)|
|[ON THE CONSTRUCTION OF  A GENERALIZED COMPUTATIONAL EXPERIMENT IN VERIFICATION PROBLEMS](https://lppm3.ru/files/journal/XLVIII/MathMontXLVIII-Alekseev.pdf):  **Article**  |![#D flow around cone](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Cone-exact-3D.png)|
|[Enhanced Pressure Based Coupled Algorithm to Combine with Pressure–Velocity-Enthalpy for all Mach Number Flow](https://link.springer.com/article/10.1007/s42405-020-00337-9):  **Article**  |![Original pressure based p–h coupled algorithm](https://media.springernature.com/lw685/springer-static/image/art%3A10.1007%2Fs42405-020-00337-9/MediaObjects/42405_2020_337_Fig1_HTML.png?as=webp)|
|[Pressure-Based Solution Framework for Non-Ideal Flows at All Mach Numbers](https://link.springer.com/chapter/10.1007/978-3-030-49626-5_4):  **Article**  |![Fully developed jet structure of a n-hexane jet injected into a quiescent nitrogen atmosphere](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-49626-5_4/MediaObjects/492738_1_En_4_Fig3_HTML.png)|
|[Entwicklung eines Simulationsmodells für Schaltlichtbögen in Überspannungsableitern](https://publikationsserver.tu-braunschweig.de/servlets/MCRFileNodeServlet/dbbs_derivate_00047899/Diss_Sander_Christian.pdf):  **PhD Thesis**  |![Spark](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/spark-tubraunscheig.png)|
|[A pressure-based solution framework for sub- and supersonic flows considering real-gas effects and phase separation under engine-relevant conditions](https://doi.org/10.1016/j.compfluid.2020.104452): **Article**|![Jet visualization](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/1-s2.0-S0045793020300281-gr15.jpg)|
|[Mixing and Autoignition of Underexpanded Methane Jets at High Pressure Conditions](https://www.research-collection.ethz.ch/handle/20.500.11850/474652): **PhD Thesis**, Results of computations with hybrid approach were used as reference for STAR-CCM|![CFD and CMC approaches](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/CFD-and-CMC.png)|

## <p align="center">  >>>>> 2019 <<<<< </p>
    
| Title | Description |
|------|-------------|
|[Numerical method to simulate detonative combustion of hydrogen-air mixture in a containment](https://doi.org/10.1080/19942060.2019.1660219):  **Article**  | ![Containement](https://www.tandfonline.com/na101/home/literatum/publisher/tandf/journals/content/tcfm20/2019/tcfm20.v013.i01/19942060.2019.1660219/20191106/images/medium/tcfm_a_1660219_f0018_oc.jpg)|
|[Numerical investigation of the auto-ignition of transient hydrogen injection in supersonic airflow](https://doi.org/10.1016/j.ijhydene.2019.07.215):  **Article**  |![Shadowgraph of the jet](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/1-s2.0-S0360319919328514-gr4.jpg)|
|[Verification on the Ensemble of Independent Numerical Solutions](https://link.springer.com/chapter/10.1007/978-3-030-22750-0_25):  **Article**  |![Solvers comparison](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-22750-0_25/MediaObjects/485772_1_En_25_Fig2_HTML.png)|
|[Computational Study of Reactants Mixing in a Rotating Detonation Combustor Using Compressible RANS](https://link.springer.com/article/10.1007/s10494-019-00097-x):  **Article**  |![Detailed shock structure of the baseline flow case in the injection region derived from the Mach number contour plot at the longitudinal mid-plane](https://media.springernature.com/lw685/springer-static/image/art%3A10.1007%2Fs10494-019-00097-x/MediaObjects/10494_2019_97_Fig5_HTML.png?as=webp)|
|[Numerical investigation of the flow characteristics of underexpanded methane jets](https://doi.org/10.1063/1.5092776): **Article**|![Methane jets flow visualization](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/methane-jet.jpeg)|
|[The numerical simulation of compressible jet at low Reynolds number using OpenFOAM](https://doi.org/10.1051/e3sconf/201912810008): **Article**|![Q-criterion for Re3600 Ma0.9 jet](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Re3600M0.9-jet.jpeg)|
|[Rocket plume URANS simulation using OpenFOAM](https://doi.org/10.1016/j.rineng.2019.100056): **Article**|![Rocket plume shadowgraph](https://ars.els-cdn.com/content/image/1-s2.0-S2590123019300568-gr3.jpg)|
|[Numerische Modellierung und Untersuchung der Hochdruckeindüsung nicht-idealer Fluide bei überkritischen Druckverhältnissen (in German)](https://athene-forschung.unibw.de/130039): **PhD Thesis**|![Jets with separation](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/overexpjets-sep.png)|

## <p align="center"> >>>>> 2018 <<<<< </p>

| Title | Description |
|------|-------------|
|[Real-Gas Effects and Phase Separation in Underexpanded Jets at Engine-Relevant Conditions](https://doi.org/10.2514/6.2018-1815):  **Article**  |![Jet development history](https://www.researchgate.net/profile/Christoph-Traxinger/publication/322309300/figure/fig5/AS:622107033612289@1525333283581/figure-fig5_W640.jpg)|
|[Analysis of the Accuracy of OpenFOAM Solvers for the Problem of Supersonic Flow Around a Cone](https://link.springer.com/chapter/10.1007/978-3-319-93713-7_18):  **Article**  |![Cone sketch](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-319-93713-7_18/MediaObjects/469704_1_En_18_Fig1_HTML.gif)|
|[Development of a new OpenFOAM solver using regularized gas dynamic equations](https://doi.org/10.1016/j.compfluid.2018.02.010):  **Article**  |![Ladenburgh jet](https://ars.els-cdn.com/content/image/1-s2.0-S0045793018300641-gr14.jpg)|
|[A hybrid pressure-based solver for nonideal single-phase fluid flows at all speeds](https://doi.org/10.1002/fld.4512):  **Article**  |![Experiment vs. calculation](https://onlinelibrary.wiley.com/cms/asset/16108f70-6fec-4197-9aab-f84cbc5c2a1d/fld4512-fig-0005-m.jpg)|
|[Comparison of the Performance of Open-Source and Commercial CFD Packages for Simulating Supersonic Compressible Jet Flows](https://doi.org/10.1109/IVMEM.2018.00019):  **Article**  |![Airbag computational domain sketch](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Airbag.png)|
|[Numerical modelling of two-dimensional perfect gas flows using RKDG method on unstructured meshes](https://doi.org/10.1063/1.5065323):  **Article** |![RKDG (a) vs. rhoPimpleCentralFoam (b)](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/forwardStep-RKDG-vs-RPCF.png) |
|[CFD methodologies for compressible atomising and cavitating multi-phase flows](https://eprints.utas.edu.au/28677/):  **PhD Thesis**  |![Jet iso contours](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Jet-Hongjiang.png)|

## <p align="center"> >>>>> 2017 <<<<< </p>

| Title | Description |
|------|-------------|
|[Numerical investigation on an array of Helmholtz resonators for the reduction of micro-pressure waves in modern and future high-speed rail tunnel systems](https://doi.org/10.1016/j.jsv.2017.04.022):  **Article** | ![Helmholtz resonantors array mesh](https://ars.els-cdn.com/content/image/1-s2.0-S0022460X17303280-gr10.jpg)|
|[Comparative Study of the Accuracy for OpenFOAM Solvers](https://doi.org/10.1109/ISPRAS.2017.00028):  **Article** |![Flow around the cone](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/index.jpeg)|
|[ Analysis of Radiation Discretization for Modelling a Spark Gap for Surge Currents ](https://doi.org/10.14311/ppt.2017.1.56):  **Article** |![Spark sketch](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Plasma-spark-sketch.png)|
|[Numerical analysis of cavitation about marine propellers using a compressible multiphase VOF fractional step method](https://www.researchgate.net/publication/319306852_Numerical_analysis_of_cavitation_about_marine_propellers_using_a_compressible_multiphase_VOF_fractional_step_method):  **Article** |![Comparison of cavitation morphology between cFSMVOF and other commercial/open-source codes](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/cFSMVOF.png)|
|[Computational analysis and mitigation of micro-pressure waves in high-speed train tunnels](https://doi.org/10.25560/72653): **PhD Thesis**|![Typical high speed train](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Train-1.png)|
|[Numerical study of characteristic modes and frequencies of flow in high-speed compressors (in English)](https://doi.org/10.15514/ISPRAS-2017-29(1)-2):  **Article** |![ERCOFTAC Pump sketch](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/ERCOFTAC-pump.png)|
|[Implementation of the solver for coupled simulation for heat transfer in gas and solid - 12th OpenFOAM Workshop](https://www.researchgate.net/publication/320924871_Implementation_of_the_solver_for_coupled_simulation_for_heat_transfer_in_gas_and_solid_-_12th_OpenFOAM_Workshop):  **Presentation** |![Problem statement](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/cht-supersonic-cone.png)|
 |[Numerical modelling of compressible flows with hybrid approximation of convective fluxes (In Russian)](https://keldysh.ru/council/3/D00202403/kraposhin_diss.pdf): **PhD Thesis**|![Comparison of density and pressure based approaches](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Dens-rho-vs-pres.png)|
|[Application of open-source software for industrial problems of vehicle lift-off gas dynamics (in Russian)](https://journals.ssau.ru/vestnik/article/view/5610): **Article**|![Comparison of Eggers experiment with calculation](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Eggers-exp-vs-calc.png)|
|[Optimization for Internal Turbulent Compressible Flows Using Adjoints](https://www.researchgate.net/publication/318144074_Optimization_for_Internal_Turbulent_Compressible_Flows_Using_Adjoints):  **Article** |![Porosity and velocity fields](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/OptIntCompFlows_AIAA.png)|

## <p align="center">  >>>>> 2016 <<<<< </p>

| Title | Description |
|------|-------------|
|[On the Stability of Supersonic Boundary Layers with Injection](https://thesis.library.caltech.edu/9755/):  **PhD Thesis** | ![Scheme of boundary layer interaction with jet](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/boundary-layer.png)|
|[Study of capabilities of hybrid scheme for advection terms approximation in mathematical models of compressible flows (in Russian)](https://ispranproceedings.elpub.ru/jour/article/view/121):  **Article** |![Liquid ring vacuum pump](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/LRVP.png)|
|[LES Discretization Methods for Unstructured Meshes Based on the Finite Volume Method](https://doi.org/10.21656/1000-0887.370228):  **Article** |![Vorticity: Flow around cylinder](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/Figs/Vorticity-vs-scheme.png)|

# For citation
[To the contents](#Contents)
    
   When using these solvers, please cite the following works:
   * [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3878441.svg)](https://doi.org/10.5281/zenodo.3878441)
   * Kraposhin MV, Banholzer M, Pfitzner M, Marchevsky IK. A hybrid pressure‐based solver for nonideal single‐phase fluid flows at all speeds. Int J Numer Meth Fluids. 2018;88:79–99. https://doi.org/10.1002/fld.4512
   * Kraposhin MV, Strijhak SV, Bovtrikova A Adaptation of Kurganov-Tadmor Numerical Scheme for Applying in Combination with the PISO Method in Numerical Simulation of Flows in a Wide Range of Mach Numbers. Procedia Computer Science. 2015;66:43-52. https://doi.org/10.1016/j.procs.2015.11.007
