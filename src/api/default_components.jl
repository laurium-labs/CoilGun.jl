module CoilGunDefaults
using CoilGun: Barrel, CoilGenerator, resistance,  Scenario, IronProjectile, ProjectilePhysical, ProjectileMagnetic, ProjectileCoilEvent, solveScenario
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn, Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Acceleration, Time, Velocity

struct NumberOfCoils
    numberOfCoils::Int
end
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
projlength = 23mm |> m
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
magIrr = 1A/m
resistor = 10Ω
volts = 15V



const default_barrel = let
    innerRadius = 1mm |> m
    thickness = 10mm |> m
    length= 0.5m |> m
    Barrel(
        innerRadius,
        thickness,
        length
        )
end

    const innerRadius = projrad+10mm 
    const outerRadius = 22.3mm |> m        #This governs how many layers of wires will be on the coil
    const length = 11.15mm |> m 
    const wireRadius = 1.6mm |> m      #This includes the insulation layer

const default_number_of_coils = let
    numberOfCoils=15
    NumberOfCoils(
        numberOfCoils
    )
end


phys    = ProjectilePhysical(projrad,
                            projlength,
                            densityFe)
mag     = ProjectileMagnetic(domainSizeFe,
                            α,
                            saturationMagnetization,
                            reversibility)
ip      = IronProjectile(phys,mag)

coils = CoilGenerator(default_number_of_coils.numberOfCoils, projrad, projrad+default_barrel.thickness, length, wireRadius)
PCE = ProjectileCoilEvent()
PCE.entersActiveZone = [nothing for _ in coils]
PCE.exitsActiveZone = [nothing for _ in coils]
totalΩ = resistor + resistance(coils[1])

endTime = 1s

default_scenario = Scenario(
    ip,
    default_barrel,
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
function solve_senario(scenario::Scenario)
    return solveScenario(scenario)
end
 end