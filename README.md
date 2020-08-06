# CoilGun.jl
## Overview ##
This program models the behavior of a ferromagnetic projectile under the influence of magnetic induction. This device, known as a Gauss Cannon (or Coil Gun), uses a series of solenoids to accelerate a projectile by sequentially turning solenoids on and off. With the aid of Julia and its packages, all known phenomena involved with magnetizing and accelerating a ferromagnetic projectile can be be calculated and simulated.

## Considerations (Accuracy of Simulation) ##
This code takes into consideration the following effects:
- Friction
- Air Resistance
- Time for switches and coils to reach a steady state value
- Energy lost in circuit and in projectile
- Primary forces as well as secondary; Secondary forces are dependent upon a dynamic system. If a system is static, no secondary forces exist. An example would be air resistance. You only feel it at high speeds.
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
## Background Theory ##
*Assuming the reader knows about basic circuital terminology, electro-magnetic physics, and knows a bit of calculus.

### Circuit/Magnetic Field
#### Current
Current is governed by the resistance in the circuit, and the voltage applied to the circuit. There are several factors in this simulation that affect the current moving through a coil at any given time. Natrually when flipping the switch to turn the current on, there will either be a steady voltage (battery) or a variable voltage (capacitor). Besides the main source of voltage, another source of voltage that is very important in consideration is the voltage induced by the moving magnetized projectile. This source of voltage opposes the main voltage source and is proportional to the projectile's velocity. What this means is that at a certain speed, the voltage induced by the projectile will equalize and maybe overcome the source voltage, effectively turning the coil off. The other factor that contributes to the voltage in the coil is self-induction and mutual induction described in the next section. Combining all of these factors into a single differential equation we obtain the following expression for voltage:
$$V = IR + L\frac{dI}{dt}+M^2\frac{d^2I}{dt^2}$$
Where $V,I,R,L,M$ are the voltage, current, resistance, self-inductance, and mutual inductance respectively. Where $M = Lk$, and $k$ is the coupling coefficient between neighboring coils. Solving for $I$ we obtain the following expression for turning the coil on:
$$I(t) = \frac{V_0}{R}(1-e^{\frac{-t}{\tau}})$$
Where $V_0$ is the supplying voltage, $t$ is time, and $\tau = \frac{Lk^2}{R}$ is the characteristic time.
For turning an operating coil off:
$$I(t) = \frac{V_0}{R}e^{\frac{-t}{\tau}}$$
#### Induction
Induction occurs when a voltage in a circuit is created by changing magnetic fields. Ultimately the Lorentz Force is to blame. This effect be caused by anything that causes a change in the magnetic field, or how the magnetic field is experienced. Some examples are a moving magnet, a moving conductive loop, other coils, or a coil acting upon itself. This effect shows up in this problem in three different ways. The first being self-inductance where the magnetic field that is being created impededes it's own creation (self-inductance), the second is an adjacent coil who's current induces a counter-current in the origional coil (mutual inductance), and lastly is the magnetized projectile that is traveling through the coils.

self-inductance is given by the relation
#### Magnetic Field
The magnetic field that is used in this simulation accounts for the length and thickness of the coil, but evaluates the magnetic field along a single axis (Z-Axis). The derivation for the equation starts with the equation for a single loop: 
$$B = \frac{\mu_0I}{2}\frac{R^2}{(R^2 + Z^2)^{\frac{3}{2}}}$$
Where $\mu_0$ is the permeablility of free-space, $I$ is the current traveling through the loop, $R$ is the radius of the loop, and $Z$ is the location along the Z-Axis with reference to the center of the loop.
From this equation, an equation that describes a coil with only a single loop can be derived by replacing $I$ with $dI$ because current is now spread out over an area. $\vec{K}$ the surface area current density which is used to describe surface current, so we will incorperate that into our equation. Here $\vec{K} \equiv \frac{d\vec{I}}{dl} = \frac{d\vec{I}}{dz} = \frac{I}{L}$, where $dl$ is the infesesimal length perpendicular to current flow. So in our case the Z-Axis is perpendicular to current flow. This gives the relation $d\vec{I} = \frac{Idz}{L}$. This in turn gives the integral:
$$\frac{\mu_0I}{2L}\int_{-a}^a\frac{dz}{(R^2+z^2)^{3/2}}$$
Giving the following equation for a solenoid with only 1 layer of windings, where $L$ is the length of the coil:
$$B(z) = \frac{\mu_0NI}{L}[\frac{z+a}{\sqrt{(z+a)^2+R^2}} - \frac{z-a}{\sqrt{(z-a)^2+R^2}}]$$
Where N is the number of windings. A similar process from getting the single loop into the single layer coil can be used to obtain a full 3-D coil. Using the current volume density $\vec{J} \equiv \frac{d\vec{I}}{da}$ to relate $\vec{I}$ and $ \vec{K}$. Where $da$ is the infintesimal area perpendicular to current flow. The relation is $\vec{J} = \frac{I}{LR} = \frac{d\vec{K}}{dr} \Rightarrow d\vec{K} = \frac{Idr}{LR}$. This gives the following integral:
$$B(z) = \frac{\mu_0NI}{LR}\int^{r_2}_{r_1}[\frac{z+a}{\sqrt{(z+a)^2+r^2}} - \frac{z-a}{\sqrt{(z-a)^2+r^2}}]dr$$
Where $r_1$ & $r_2$ were the inner and outer radius of the coil and $R$ is now $r_2 - r_1$ or the thickness of the coil. This integration gives the following equation for a coil of uniform volume current density:
$$B(z) = \frac{\mu_0NI}{LR}[(z+a)ln(\frac{\sqrt{(z+a)^2+r_2^2}+r_2}{\sqrt{(z+a)^2+r_1^2}+r_1})-(z-a)ln(\frac{\sqrt{(z-a)^2+r_2^2}+r_2}{\sqrt{(z-a)^2+r_1^2}+r_1})]$$

### Magnetization
The magnetization of a material indicates how aligned it's magnetic dipoles are in a single direction. For most materials, the magnetization is a linear relationship between the applied field, and the resulting magnetization. For ferromagnetic materials, the relationship is a bit more complex. The magnetization of a ferromagnetic material depends upon the previous magnetization of the material. This is known as hysteresis. 

Magnetization is the result a few competing factors. Namely, diamagnetism, paramagnetism, and ferromagnetism. The changing orbital path of the electron contributes to diamagnetism, the spin of the electrons and protons in the atom contributes to paramagnetism, and the group alignments of the magnetic dipoles of a group of atoms contributes to ferromagnetism.

When a ferromagnetic object is unmagnetized, there exists small groups of atoms who's magnetic dipoles are all aligned in a random direction called magnetic domains. These domains are all roughly the same size, and share a boundry of atoms whose dipoles adjust from one domian direction to another called a domain wall. These walls vary in size, but are much smaller than the domian size itself. The reason why a ferromagnetic material can be unmagnetized is due to these domians pointing in randoms directions, effectively canceling each other's magnetization out. When a ferromagnetic object undergoes magnetization for the first time the domains that are aligned with the applied magnetic field will expand. This expansion is caused the movement of the domain walls. 

A macroscopic way to represent domain wall movement is gather a long line of people. The first person will reach upward and leave the hands upward. The second person will follow in suit, copying the first person as fast as possible. Likewise the third person will copy the second and so on. The resulting "wave" of hands is representing the reorientation of individual atomic magnetic dipoles changing to the new direction, and the amount of hands raised is the growing magnetic domian.
#### Reversible Magnetization
Magnetization can be thought of as two different processes that happen inside the material, reversible magnetization and irriversible magnetization. As the name suggest reversible magnetization is a type of magnetization that can be reversible. What this means is that an object returns back to an unmagnetized state after a magnetic field is appied and then removed. This can happen for ferromagnetic materials under small (compared to it's saturation magnetization) magnetic fields, or materials that have no or very little crystal defects. Crystal is used in the sense of atoms arrayed periodically forming a lattice.

This type of magnetism can be represented with linear relationship between the anhysteresis magnetization $M_{an}$ and the irriversible magnetization $M_{irr}$ Ref.[2]:
$$M_{rev} = c(M_{an} - M_{irr})$$
 $M_{an}$ is also called the Langevin equation described as:
$$\mathscr{L}(H_e) = M_s[coth(\frac{H_e}{a}) - \frac{a}{H_e}]$$
Where $H_e$ is the effective magnetic field experienced by each domain, here represented as: $H_e = H + \alpha M$. Where $\alpha$ is the magnetic interdomain coupling factor, typically determined from experimentation. $a$ is a constant used to normalize $H_e$ and determines the density distribution of domain clusters, and $M_s$ is the saturation magnetization. In the next section $M_{irr}$ will be described.
#### Irriversible Magnetization
Irriversible magnetization on the other hand does not return to it's origional unmagnetized state, making it a bit harder to represent mathematically. Instead there is some remaining magnetism to the object after the applied field is removed, called remenance. This remenance are caused by defects in the crystal lattice, and are reffered to as pinning sites. This is because they typically hard to alter their magnetic domain orientation than normal. This, in turn, is affecting neighboring domains as well. This remenance is caused by when domain walls encounter and overcome these pinning sites in the crystal, changing the direction of the defects dipole orientation. This then results in a total non-randomized dipole direction of the overall object. A simple equation that is used to calculate the irreversible magnetization Ref.[2]:
$$M_{irr} = \frac{M - cM_{an}}{1-c}$$
Where $M$ is the total magnetization of the material, and $c$ is the reversibility parameter. The reversiblity parameter describes how strong a magnetic field needs to be to demagnetize a magnetized object. A high $c$ value means that the object is easily magnetized.
#### Calculation of Magnetization
The magnetization of an object is the sum of both the reversible magnetization and the irriversible magnetization. Shown as $M = M_{irr} + M_{rev}$. Ref.[2]

In order to calculate the magnetizm of the object iteratively, the change in magnetism per change in the applied magnetic field needs to be calculated. This turns the above equation into:
$$\frac{dM}{dH} = \frac{dM_{irr}}{dH} + c[\frac{dM_{an}}{dH} - \frac{dM_{irr}}{dH}]$$
Because both $M_{irr}$ & $M_{an}$ are calculated using $H_e$, the components need to be changed to include $H_e$. The following relations are done as such:
$$\frac{dM_{irr}}{dH} = \frac{dM_{irr}}{dH_e}\frac{dH_e}{dH} = \frac{\delta_M(M_{an}-M_{irr})}{k\delta}[1 + \alpha \frac{dM}{dH}]$$
$$\frac{dM_{an}}{dH} = \frac{dM_{an}}{dH_e}\frac{dH_e}{dH} = \frac{d\mathscr{L}}{dH_e}[1 + \alpha \frac{dM}{dH}]$$
Where $\delta = sign[\frac{dH}{dt}]$ (meaning, if H is increasing $\delta$ is 1 and if decreasing -1). $\delta_M = \frac{1 + sign[M_{rev}\delta]}{2}$, which takes the values of 1 and 0, and $k$ is the domian pinning factor. The reason why both $\delta$ & $\delta_M$ are used is due to pinning effects. These pinning sites always resist magnetization, so shape of the curve changes as a result of the changing magnetic field, which is the cause of hysteresis for ferromagnetism. $k$, the domain pinning factor, in this simulation is treated as a constant, but it does vary. 
Summing the above two equations and solving for $\frac{dM}{dH}$ gives the following equation that is used to calculate the change in magnetism.
Due to pinning effects, the change in the magnetization curve changes depending on whether the applied magnetic field is increasing or decreasing. 
$$\frac{dM}{dH} = \frac{\frac{dM_{irr}}{dH_e} - c\frac{dM_{an}}{dH_e}}{1 - \alpha(\frac{dM_{irr}}{dH_e} - c\frac{dM_{an}}{dH_e})}$$

In order to just obtain $dM$, the amount that the magnetization changes, the above function needs to be mulitplied by $dH$. Here $dH$ is calulated by the following method:
$$dH = (\nabla H*v + \frac{\partial H}{\partial I}*\frac{\partial I}{\partial t})dt$$
Where $v$ is the velocity of the projectile, and $I$ is the current going through the coil. The reason why $dH$ is broken up into two pieces is because $dH$ is a function of both position and current.

### Dipole-Coil Interaction
The force that exists between a coil and a magnet is described in this simulation as:
$$F = \vec{m} • \nabla \vec{B}$$
What this means is that the strength of the magnetic field doesn't affect the force between the two, but how much the magnetic field changes (per unit length). $\vec{m}$ here is the magnetic moment of the projectile found through the relation:
$$\vec{m} = \int\vec{M}dV$$

### Resistive Forces
In this simulation there exists two forces that resist motion, drag and friction.
#### Drag
The method used to find the drage equation is fairly straightforward and is commonly used elsewhere. The equation used to calculate drag is:
$$F_{drag} = 6\pi\mathscr{R}rv$$
Where $\mathscr{R}$ is the dynamic viscosity of air, $r$ is the radius of the projectile, and $v$ is the velocity of the projecitle.
#### Friction
Two different types of friction are used in this simulation, static friction and dynamic friction. Static friction is applied when the object is not moving, and dynamic friction is used when the object is moving. The general friction equation is :
$$F_{friction} = \mu F_N = \mu mg$$
Where $\mu$ is the coefficient of fricition (either static of kinetic), $m$ is the mass of the projectile, and $g$ is the gravitational acceleration that the projectile experiences.


## Construction of the Coil Gun ##
### The Simulation
The most critical aspect of this simulator is the simulation of a coil. So naturally, the simulation started by creating functions on how to calculate the physical aspects of the coil. Later functions like calculating the magnetic field produced from the coil were developed along with functions relating to how the magnetization of the projectile changed with the changing magnetic fields. From there the functions that decribed how a moving magnetized projectile interacted with the coil were developed including the induced voltage, and the force between the two (along with the resistance forces). The functions that described self-inductance and mutual inductances were then created and implemented. Up until this point, there was not enough functions to simulate the general characteristics of how the the projectile could be pulled through the coil, so it wasn't until after this point where the solver/simulator was created and developed (by Brent). After a simulation of a single coil was complete, then functions and changes to existing functions were performed in order to allow for multiple coils. Once any number of coils could be used in the simulation, the functions used to calculate the force and related aspects of the simulation got increasingly more complex to better match reality.

## Applications ## 
The device itself can be used for a number of different applications. Some notable applications are: general construction, clearing computer memory, material testing, scientific demonstrations. The goal of this project is to accelerate the projectile to mach 1.
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