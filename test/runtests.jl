using CoilGun: Barrel, CoilGenerator, resistance,  Scenario, IronProjectile, ProjectilePhysical, ProjectileMagnetic, ProjectileCoilEvent, solveScenario, Projectile
using CoilGun: Coil, Barrel, ProjectileCoilEvent, Voltage, ElectricalResistance, HField, Length, Velocity
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn
using Unitful:Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Acceleration, Time, Velocity
using ForwardDiff
#using Plots
using Test
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
#Projectile Specifications
#Physical
projrad = 3.5mm |> m
projlength = 1.0inch |> m
position = 0m   
velocity = 0.1m/s
accel = 0m/s^2
#Magnetic
saturationMagnetization = 1.61e6A/m
reversibility = 0.373
const domainPinningFactor = 742.64A/m           # This is the domain pinning factor from Ref.[5]
const α = 1.34e-3                               # Interdomain Coupling Factor from Ref.[5]
const a = 882.55A/m                             # "Determines the density distribution of mag. domians"~Ref.[2] Ref.[5]
const magMomentPerDomain = k*roomTemp/a         # This dipole magnetic moment from Ref.[5]
magnetization = 0A/m
magIrr = 1.0A/m
println("Initial Parameters:\n\tposition:\t\t",position,
    "\n\tvelocity:\t\t", velocity,
    "\n\tacceleration:\t\t", accel, 
    "\n\tmagnetization:\t\t", magnetization)

#Coil Specifications
bthickness = 1mm |> m #Barrel thickness

innerRadius = projrad+bthickness        #The Inner diameter of the Coil needs to be the Outer diameter of the barrel
cthickness = 1inch |> m                 #The difference in the Inner diameter and Outer diameter of the Coil
coilLen = 1.0inch |> m                  #The length of the Coil should be the exact length of the projectile
coilHght = 2.3e-2m |> m                 #Distance from inner to outer diameter of the Coil
wirerad = 1.6mm |> m                    #The radius of 14-guage wire including insulation
resistor = 10Ω
volts = 15V
numberOfCoils = 4

#Barrel Specifications
blength = numberOfCoils * coilLen |> m #Length of barrel
bRadius = projrad

phys    = ProjectilePhysical(projrad,
                            projlength,
                            densityFe)
mag     = ProjectileMagnetic(domainSizeFe,
                            α,
                            saturationMagnetization,
                            reversibility)
ip      = IronProjectile(phys,mag)
bar     = Barrel(bRadius,bthickness,blength)
coils = CoilGenerator(numberOfCoils, innerRadius, innerRadius+cthickness, coilLen, wirerad)
PCE = ProjectileCoilEvent()
PCE.entersActiveZone = [nothing for _ in coils]
PCE.exitsActiveZone = [nothing for _ in coils]


# Δt=0.1s
# t = 0s
# Magirr = 0A/m
# coil = coils[1]
totalΩ = resistor + resistance(coils[1])
# Curr = map(i ->CoilGun.current(coils[i], totalΩ, volts, t - eventTimes[1].entersActiveZone, magnetization, velocity, position), 1:length(coils))
# I = Curr[1]
# B = bFieldCoil(coil, I, position)
# ∇B = ∇BFieldCoil(coil, I, position)
# Magirr = Mag_irr(ip, B, Magirr, magnetization)
# dH = dHField(coils, volts, totalΩ,∇B, magnetization, position, velocity, accel, t) 
# magnetization += ∂Magnetization_∂HField(ip, B, Magirr, dH) * dH * Δt

# ∂current = ∂Current(coil, t, volts, totalΩ, position, velocity, acceleration(totalForce(ip, ∇B, velocity, magnetization), mass(ip)), magnetization)
# ∂BField_∂Current(coil, position)
# include("unitTests.jl")
# return
# println("Finished Unit tests")

endTime = 0.2s

scenario = Scenario(
    ip,
    bar,
    coils,
    PCE,
    endTime,
    volts,
    resistor,
    magIrr,
    position,
    velocity,
    magnetization
)

@show scenario
println("Now Solving...")
sln = solveScenario(scenario)
println("Length of Velocity:\t\t",length(sln[4,:]),"\nLength of Position:\t\t", length(sln[3,:]),"\nLength of Time:\t\t\t", length(sln[2,:]),"\nLength of Magnetization:\t", length(sln[1,:]))
xAxis = 1:length(sln[1,:])
#figure(figsize=(8,6))
p1 = plot(sln, vars=(0,2), title = "Displacement", ylabel = "[m]")
p2 = plot(sln, vars=(0,3), title = "Velocity", ylabel = "[m/s]")
p3 = plot(sln, vars=(0,1), title = "Magnetization", ylabel = "[A/m]", legend = false)
p4 = plot(sln, vars=(0,4), title = "Irriversible Magnetization", ylabel = "[A/m]")
display(plot(p1,p2,p3,p4, layout = (2,2)))

# dist = coilLen|>m|>ustrip
# println("Max Velocity $(sln[3,:][argmax(sln[3,:])])m/s")
# println("Point were projectile started to accerate $(dist .- sln[2,:][argmin(sln[3,:])])m.\nPoint where projectile started to decelerate $(sln[2,:][argmax(sln[3,:])] .- dist)m")
# plot(sln, layout = (2,2))