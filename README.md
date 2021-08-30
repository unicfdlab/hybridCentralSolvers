# hybridCentralSolvers
United collection of hybrid  Central solvers - one-phase, two-phase and multicomponent versions.

Only OpenFOAM+ version of the OpenFOAM technology is supported sinces 2018. Use branches digitef-dev-YYMM, where YYMM - corresponds to OpenFOAM+ version, for example 1812

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

Telegram public group: https://t.me/hybridCentralSolvers for questions and discussion of the solvers features

Available OpenFOAM versions:
* OpenFOAM 3.1 - master branch
* OpenFOAM 4.1 - dev-of4.1 branch
* OpenFOAM 6   - dev-of6 branch
* OpenFOAM+ 1812 - digitef-dev-1812
* OpenFOAM+ 1912 - digitef-dev-1912
* OpenFOAM+ 2012 - digitef-dev-2012

## The library was useful in next research studies.

### >>>>> 2021 <<<<<

| Title | Description |
|------|-------------|
|[Development of a new OpenFOAM solver using regularized gas dynamic equations](https://www.mdpi.com/2311-5521/6/8/274)| ![Jet with particles Logo](https://www.mdpi.com/fluids/fluids-06-00274/article_deploy/html/images/fluids-06-00274-ag-550.jpg) |

###  >>>>> 2019 <<<<<
| Title | Description |
|------|-------------|
|[Numerical method to simulate detonative combustion of hydrogen-air mixture in a containment](https://doi.org/10.1080/19942060.2019.1660219)| ![Containement](https://www.tandfonline.com/na101/home/literatum/publisher/tandf/journals/content/tcfm20/2019/tcfm20.v013.i01/19942060.2019.1660219/20191106/images/medium/tcfm_a_1660219_f0018_oc.jpg)|
|[Numerical investigation of the auto-ignition of transient hydrogen injection in supersonic airflow](https://doi.org/10.1016/j.ijhydene.2019.07.215)|![Shadowgraph of the jet](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/1-s2.0-S0360319919328514-gr4.jpg)|
|[Verification on the Ensemble of Independent Numerical Solutions](https://link.springer.com/chapter/10.1007/978-3-030-22750-0_25)|![Solvers comparison](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-22750-0_25/MediaObjects/485772_1_En_25_Fig2_HTML.png)|

### >>>>> 2018 <<<<<
| Title | Description |
|------|-------------|
|[Real-Gas Effects and Phase Separation in Underexpanded Jets at Engine-Relevant Conditions](https://doi.org/10.2514/6.2018-1815)|![Jet development history](https://www.researchgate.net/profile/Christoph-Traxinger/publication/322309300/figure/fig5/AS:622107033612289@1525333283581/figure-fig5_W640.jpg)|
|[Analysis of the Accuracy of OpenFOAM Solvers for the Problem of Supersonic Flow Around a Cone](https://link.springer.com/chapter/10.1007/978-3-319-93713-7_18)|![Cone sketch](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-319-93713-7_18/MediaObjects/469704_1_En_18_Fig1_HTML.gif)|

### >>>>> 2017 <<<<<
| Title | Description |
|------|-------------|
|[Numerical investigation on an array of Helmholtz resonators for the reduction of micro-pressure waves in modern and future high-speed rail tunnel systems](https://doi.org/10.1016/j.jsv.2017.04.022)| ![Helmholtz resonantors array mesh](https://ars.els-cdn.com/content/image/1-s2.0-S0022460X17303280-gr10.jpg)|

###  >>>>> 2016 <<<<<
| Title | Description |
|------|-------------|
|[On the Stability of Supersonic Boundary Layers with Injection](https://thesis.library.caltech.edu/9755/)| ![Scheme of boundary layer interaction with jet](https://github.com/unicfdlab/hybridCentralSolvers/blob/master/boundary-layer.png)|

   When using these solvers, please cite the following works:
   * [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3878441.svg)](https://doi.org/10.5281/zenodo.3878441)
   * Kraposhin MV, Banholzer M, Pfitzner M, Marchevsky IK. A hybrid pressure‐based solver for nonideal single‐phase fluid flows at all speeds. Int J Numer Meth Fluids. 2018;88:79–99. https://doi.org/10.1002/fld.4512
   * Kraposhin MV, Strijhak SV, Bovtrikova A Adaptation of Kurganov-Tadmor Numerical Scheme for Applying in Combination with the PISO Method in Numerical Simulation of Flows in a Wide Range of Mach Numbers. Procedia Computer Science. 2015;66:43-52. https://doi.org/10.1016/j.procs.2015.11.007
