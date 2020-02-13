# CoilGun.jl
This program models the acceleration of an iron prjectile with the use of magnetic induction. This device, known as a Gauss Cannon (or Coil Gun), uses a series of electromagnetic coils for acceleration. This code takes into consideration the following effects: Friction, air resistance, changing magnitic domain sizes (100X oversized domain sizes), time taken to turn on switches and coil, energy lost in circuit and in projectile, and finally primary forces as well as secondary. Secondary forces are dependent upon a dynamic system. If a system is static, no secondary forces exist. An example would be the resistance to current when turning a electromagnetic coil on.

This project assumes* that the magnetic field inside a wire loop is constant, and is calculated at the center. It's not taking any quantum effects (besides the domain magnetic dipoles) into consideration. Not taking the frictional force of rifling into consideration or the effects that it would have on the projectile, only assuming Newtonian drag and not Navier-Stokes, not shure what will be used for air resistance in the barrel, might be taking the capacitance of the coil into consideration, not taking earths magnetic field into consideration, and the projectiles geometry is a perfect cyliner.

*All assumptions made will gradually be made to fit closer with reality. The complexity of the system will increase until either it at at the highest complexity, and matches reality, or until computation time is noticably effected (like taking 1 minute to run). This is the reason why the magnetic domain size of the projectile will not be able to fit closer with reality, and thus will be at this level of approximation.





```math
```