# hybridCentralSolvers
United collection of hybrid  Central solvers - one-phase, two-phase and multicomponent versions.

Latest OpenFOAM version: 4.1

1. **pimpleCentralFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
Kurganov and Tadmor with LTS support for steady-state calculations

2. **pimpleCentralDyMFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
Kurganov and Tadmor with mesh motion and LTS support for steady-state calculations

3. **reactingPimpleCentralFoam** - Pressure-based semi implicit compressible flow solver based on central-upwind schemes of 
Kurganov and Tadmor for combustion with chemical reactions and LTS support for steady-state calculations

4. **twoPhaseMixingCentralFoam** - Transient Eulerian two-phase solver. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface is not accounted.

5. **twoPhaseMixingCentralDyMFoam** - Transient Eulerian two-phase solver with dynamic meshes. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface is not accounted.

6. **chtMultiRegionCentralFoam** -     Pressure-based semi implicit solver, based on hybrid central-upwind schemes
    of Kurganov and Tadmor for conjugate simulation of compressible flows (Mach 
    number is ranging from 0 to 6) and solid body heat transfer.