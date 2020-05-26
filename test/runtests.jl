using CoilGun
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn
using Unitful:Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, 
                Inductance, ustrip, Voltage, Acceleration, Time, Velocity
using ForwardDiff
using Plots

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
const roomTemp = 293K                                       #Standard room Tempearture
const domainMagnetization = 0.2 * numberAtomsperDomainFe*bohrMagnetonPerAtomFe |> A/m #Magnetization of the domain
const saturationMagnetizationPerKgFe = 217.6A/(m*kg)             #Saturation magnetizaiton of pure Iron per unit mass.
const kineticFrictionCoefficientFe = 0.36                   #Kinetic friction coefficient of Mild Steel on Copper, probably not exact
const staticFrictionCoefficientFe = 0.53                    #Static friction coefficient of copper on Steel, probably not exact
const dynamicViscosityAir = 1.825e-5kg/(m*s)                #Dynamic viscosity of air at 20C

magnetization = domainMagnetization
#Projectile Specifications
#Physical
projrad = 3.5mm |> m
projlength = 1.0inch |> m
position = projlength/2
velocity = -0.1m/s
accel = 0m/s^2
#Magnetic
saturationMagnetization = 1.61e6A/m
reversibility = 0.373
const domainPinningFactor = 742.64A/m           #This is the domain pinning factor from Ref.[5]
const α = 1.34e-3                               #Interdomain Coupling Factor from Ref.[5]
const a = 882.55A/m                             #"Determines the density distribution of mag. domians"~Ref.[2] Ref.[5]
const magMomentPerDomain = k*roomTemp/(a * μ0)  #This dipole magnetic moment from Ref.[5]


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
totalΩ = resistor + resistance(coil)
I = CoilGun.current(ip, coil, totalΩ, volts, t, magnetization, velocity, position)
B = bFieldCoil(coil, I, position)
∇B = ∇BFieldCoil(coil, I, position)
Magirr = Mag_irr(ip, B, Magirr, magnetization)
dH = ∂HField(coil, I, volts, totalΩ,∇B, magnetization, position, velocity, accel, t) 
magnetization += ∂Magnetization_∂HField(ip, B, Magirr, dH) * dH * Δt

∂current = ∂Current(coil, t, volts, totalΩ, position, velocity, acceleration(totalForce(ip, ∇B, velocity, magnetization), mass(ip)), magnetization)
∂BField_∂Current(coil, I, position)

endTime = 0.2s

scenario = Scenario(
    ip,
    bar,
    coil,
    endTime,
    volts,
    resistor,
    Magirr,
    position,
    velocity,
    magnetization
)
sln = solveScenario(scenario)

p1 = plot(sln, vars=(0,2), title = "Position")
p2 = plot(sln, vars=(0,3), title  = "Velocity")
p3 = plot(sln, vars=(0,1), title  = "Magnetization")
p4 = plot(sln, vars=(0,4), title  = "Irriversible Magnetization")
plot(p1,p2,p3,p4)