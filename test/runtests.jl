using CoilGun
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn
using Unitful:Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, 
                Inductance, ustrip, Voltage, Acceleration, Time, Velocity
using ForwardDiff

const resistivityCu = 1.72e-8m*Ω                            #Resistivity of Copper
const densityCu = 8960kg/m^3                                #Density if pure Copper
const densityFe = 7750kg/m^3                                #Density of pure Iron
const atomicWeightFe = 55.845g/mol  |> kg/mol               #Atomic weight of Iron
const domainSizeFe = 26.5e-7m                               #Average magnetic domain size for pure Iron (actually 26.5e-9)
const densityNi = 8.908g/cm^3 |> kg/m^3                     #Density of pure Nickel
const currieTempFe = 1043K                                  #This is the Currie tempearture of Iron
const bohrMagnetonPerAtomFe = 2.2m^-3*μB   |> A/m           #The is the Bohr Magneton per Iron atom
const numberAtomsperDomainFe = Na*domainSizeFe^3*densityFe/atomicWeightFe  #Number of atoms per Iron domain volume
const magPerFeAtom = currieTempFe*k/bohrMagnetonPerAtomFe   #This is the magnetic field given off by each Iron atom
const magPerFeDomain = magPerFeAtom*numberAtomsperDomainFe  #Magnetic field of the domain
const χFe = 200_000                                         #Magnetic susceptibility of iron at 20 C (unitless)
const μ = μ0*(1+χFe)                                        #Magnetic pearmeability of iron
const α = 9.5e-5                                            #Interdomain Coupling Factor (for an iron transformer)
const roomTemp = 293K                                       #Standard room Tempearture
const domainPinningFactor = 150A/m                          #This is the domain pinning factor for Iron (transformer)
const domainMagnetization = 0.2 * numberAtomsperDomainFe*bohrMagnetonPerAtomFe |> A/m #Magnetization of the domain
const magMomentPerDomain = domainMagnetization*domainSizeFe^3    #This dipole magnetic moment doesn't take hysteresis/pinning into effect
const saturationMagnetizationPerKgFe = 217.6A/(m*kg)             #Saturation magnetizaiton of pure Iron per unit mass.

magnetization = domainMagnetization
#Projectile Specifications
projrad = 3.5mm |> m
projlength = 1inch |> m
magdomiansize = 26.5μm |> m   #From literature, it is actually 26.5 nm, this value is not being used because I would have no memory left.
magstrngth = 1T
saturationMagnetization = saturationMagnetizationPerKgFe * magdomiansize^3 * densityFe
position = projlength
velocity = 1m/s
accel = 0m/s^2
reversibility = 0.373

#Barrel Specifications
bthickness = 1mm |> m #Barrel thickness
blength = 0.5m |> m #Length of barrel


#Coil Specifications
innerRadius = projrad+bthickness        #The Inner diameter of the Coil needs to be the Outer diameter of the barrel
cthickness = 1inch |> m                 #The difference in the Inner diameter and Outer diameter of the Coil
coilLen = projlength                    #The length of the Coil should be the exact length of the projectile
coilHght = 2.3e-2m |> m                 #Distance from inner to outer diameter of the Coil
wirerad = 1.6mm |> m                    #The radius of 14-guage wire including insulation
I = 1A                                  #Current flowing through the wire
stepSize = 1_000
resistor = 10Ω
volts = 15V

phys    = ProjectilePhysical(projrad,
                            projlength,
                            densityFe)
mag     = ProjectileMagnetic(domainSizeFe,
                            α,
                            saturationMagnetization,
                            reversibility)
ip      = IronProjectile(phys,mag)
bar     = Barrel(ip.physical.radius,bthickness,blength)
coil    = Coil(projrad,projrad+cthickness,ip.physical.length,wirerad)


Δt=0.1s
t = 2s
Magirr = 0A/m
prevI = 0A
totalΩ = resistor + resistance(coil)
I = current(ip, coil, totalΩ, volts, t, magnetization, velocity, position)

B = simpleBField(coil, I, position)
∇B = bFieldGradient(coil, I, position)
Magirr = Mag_irr(ip, B, Magirr, magnetization)
dH = ∂HField(coil, I, volts, totalΩ,∇B, magnetization, position, velocity, accel, t) 
# (∇B * ip.velocity + simpleBField(coil, I - prevI, ip.position)/t) * t / μ0
magnetization += ∂Magnetization_∂HField(ip, B, Magirr, dH) * dH * Δt

∂current = ∂Current(coil, t, volts, totalΩ, position, velocity, acceleration(totalForce(ip, ∇B, velocity, magnetization), mass(ip)), magnetization)
∂SimpleBField_∂Current(coil, I, position)
