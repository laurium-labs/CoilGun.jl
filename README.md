# CoilGun.jl

## Overview ##
This program models the behavior of a ferromagnetic projectile under the influence of magnetic induction. This device, known as a Gauss Cannon (or Coil Gun), uses a series of solenoids to accelerate a projectile by sequentially turning solenoids on and off. With the aid of Julia and its packages, all known phenomena involved with magnetizing and accelerating a ferromagnetic projectile can be be calculated and simulated.

The device itself can be used for a number of different applications. Some notable applications are: general construction, clearing computer memory, material testing, scientific demonstrations, etc, etc. The goal of this project is to accelerate the projectile to mach 1.


This code takes into consideration the following effects:
- Friction
- Air Resistance
- Time for switches and coils to reach a steady state value
- Energy lost in circuit and in projectile
- Primary forces as well as secondary; Secondary forces are dependent upon a dynamic system. If a system is static, no secondary forces exist. An example would be an induced counter-current when a magnet is traveling towards or away from a solenoid.
- Pinning effects of the given material

Assumptions*:
- The magnetic field inside a wire loop is constant, and is calculated at the center. Also, all of the winding are treated as if they're in a single loop at the effective radius.
- Projectiles geometry is a perfect cyliner
- The pinning factor (when calculating the magnetizaiton) is a constant i.e. energy loss due to domain rotation is negligible.    
- The interaction between the coil and the projectile is like a wire loop to a magnetic dipole
- Each coil/loop is seperated by the actual length of the coil
- The current created by the capacitor and coupled coils is completely independent from the current induced by the accelerating projectile.
- The entire material's magnetizaiton can be assumed to be constant throughout the projectile, holding the magnetization value of the closest part of the projectile to the magnetic field. Also that the behavior of the magnetic domains can be accurately generalized and averaged.
- The friciton force between the projectile and the barrel is, once moving, independent upon velocity.


It's not taking into consideration:
- Quantum effects (besides the domain magnetic dipoles) into consideration
- Frictional force of rifling into consideration or the effects that it would have on the projectile
- Navier-Stokes drag (using Newtonian)
- Magnetostriction, other stress/strain relations to magnetizm or their effect on permeability.
- Earths magnetic field  
- How the barrel affects the generated magentic fields.
- How the voltage produced from the capacitor changes in time.
 

*All assumptions made will gradually be made to fit closer with reality. The complexity of the system will increase until either the simulation provieds accurate resutls and matches reality (error < 1%), or until computation time is unreasonably long (like taking 10 minutes to complete the simulation).

#### Sources: ####

[1] Jiles, D. C., & Atherton, D. L. (1986). Theory of ferromagnetic hysteresis. Journal of magnetism and magnetic materials, 61(1-2), 48-60.
    
   - This paper was one of the origional papers that created the currently used model. This specific paper was used to understand the magnetization process. This describes the derivation of some of the equations used. Some of the equations used were:
    Langevin equation: 
    $\mathscr{L}(\frac{Be}{a}) = coth(\frac{Be}{a}) - \frac{a}{Be}$
    Where $Be = \mu_0(H + \alpha M)$, $\alpha$ representing interdomain coupling between magnetic domains inside the projectile and $M$ being the magnetization of the entire projectile. $a = \frac{k_BT}{m}$ with $k_B$ the Boltzmann constant, $T$ is temperature in kelvin, and $m$ is the magnetization moment of a domain. The Langevin equation represents how a ferromagentic material, without defects, would respond to a magnetic field.

[2] Chwastek, K., & Szczygłowski, J. (2008). Estimation methods for the Jiles-Atherton model parameters–a review. Przegląd Elektrotechniczny, 84(12), 145-147.
- This paper provides an excellent review of the Jiles-Atherton model, and gives a method for obtaining the necessary parameters.

[3] Nicholas, R.J. (2007). Magnetic Properties of Materials. http://www-rjn.physics.ox.ac.uk/lectures/magnetismnotes10.pdf
- Used to better understand how the magnetization process works.

[4] Jarrett Revels, Miles Lubin, and Theodore Papamarkou. Forward-Mode Automatic Differentiation in Julia. 2016. arXiv:1607.07892 [cs.MS].
- This paper describes how ForwardDiff works. A package that was used in the simulation to provide very quick derivatives.

[5] Sedira, D., Gabi, Y., Kedous-Lebouc, A., Jacob, K., Wolter, B., & Straß, B. (2020). ABC method for hysteresis model parameters identification. Journal of Magnetism and Magnetic Materials, 166724.
- Used parameters given in paper for comparision

[6] Li, L. (2004). Stress effects on ferromagnetic materials: Investigation of stainless steel and nickel.
- Used for additional understanding of magnetization and implementation of a different form than Ref.[1] that was used to describe the magnetization.

[7] Zhu, B. (2001). Non-linear, irreversible magnetization processes in magnetic materials: instrumentation, measurements, modeling and application.
- Similar to Ref.[6], but the Jiles-Atherton model that was derived in this thesis was the model used in the simulation.

[8] Bastos, J., Sadowski, N. (2003). Electromagnetic Modeling by Finite Element Methods. Boca Raton: CRC Press, https://services.lib.mtu.edu:5021/10.1201/9780203911174
- Used to obtain the Taylor series approximation of the $\mathscr{L}$ function, and a greater understanding of the derivation of the $\frac{dM}{dH}$ function.