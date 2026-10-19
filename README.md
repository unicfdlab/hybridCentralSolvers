# hybridCentralSolvers
United collection of hybrid  Central solvers - one-phase, two-phase and multicomponent versions. Full list of solvers is available at the master branch.

Telegram public group: https://t.me/hybridCentralSolvers for questions and discussion of the solvers features

## Building

This branch is laid out as an OpenFOAM (openfoam.com) plugin.

- Standalone build into user locations (`$FOAM_USER_APPBIN`, `$FOAM_USER_LIBBIN`):
  ```bash
  ./Allwmake
  ```
- As a plugin: place (or add as a git submodule) this repository in
  `$WM_PROJECT_DIR/plugins/hybridCentralSolvers` and run `./Allwmake-plugins` from
  `$WM_PROJECT_DIR`. It is then installed into `$FOAM_APPBIN`, `$FOAM_LIBBIN`.
- Custom location:
  `./Allwmake -prefix=user|group|openfoam|/some/path`

   When using these solvers, please cite the following works:
   * Kraposhin MV, Banholzer M, Pfitzner M, Marchevsky IK. A hybrid pressure‐based solver for nonideal single‐phase fluid flows at all speeds. Int J Numer Meth Fluids. 2018;88:79–99. https://doi.org/10.1002/fld.4512
   * Kraposhin MV, Strijhak SV, Bovtrikova A Adaptation of Kurganov-Tadmor Numerical Scheme for Applying in Combination with the PISO Method in Numerical Simulation of Flows in a Wide Range of Mach Numbers. Procedia Computer Science. 2015;66:43-52. https://doi.org/10.1016/j.procs.2015.11.007
