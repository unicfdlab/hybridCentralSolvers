ToDo:
- add viscousity into volume fraction evolution equation
- interface (vF1) Courant
- revise the code
- check conservation of the liquid volume fraction


InProgress:
- run shockwave / droplet interaction+


Done:
- General integration scheme (algortihm)
- Simple Riemann test cases
- Dam break and gravity
- Energy equation
- surface tension force+
- run cylinder collapse+
- turbulence model+
- remove diffusion terms from phase energy equations+
- viscous dissipation in energy equation+
- U = (rhoU)_own / rho_own, where rho_own vF1*rho1_own + vF2*rho2_own,
(rhoU)_own = vF1*(rho1U)_own + vF2*(rho2U)_own, etc+ (doesn't improve the
situation too much, at least for the first attempt)
