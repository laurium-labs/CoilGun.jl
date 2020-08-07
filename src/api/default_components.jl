module CoilGunDefaults
using CoilGun: Barrel, CoilGenerator, resistance,  Scenario, IronProjectile, ProjectilePhysical, ProjectileMagnetic, ProjectileCoilEvent, solveScenario, Projectile
using CoilGun: Coil, Barrel, ProjectileCoilEvent, Voltage, ElectricalResistance, HField, Length, Velocity
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn, Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Acceleration, Time, Velocity

struct UIScenario
    ip::Projectile
    barrel::Barrel
    coil::Coil
    endTime::Time
    voltage::Voltage
    resistor::ElectricalResistance
    initialMagIRR::HField
    initalPosition::Length
    initialVelocity::Velocity
    initialMagnetization::HField
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
projlength = 2.54cm |> m
position = 0.0m   
velocity = 0.1m/s
accel = 0.0m/s^2
#Magnetic
saturationMagnetization = 1.61e6A/m
reversibility = 0.373
const domainPinningFactor = 742.64A/m           # This is the domain pinning factor from Ref.[5]
const α = 1.34e-3                               # Interdomain Coupling Factor from Ref.[5]
const a = 882.55A/m                             # "Determines the density distribution of mag. domians"~Ref.[2] Ref.[5]
const magMomentPerDomain = k*roomTemp/a         # This dipole magnetic moment from Ref.[5]
magnetization = 0.0A/m
initialMagIRR= 1.0A/m
resistor = 10Ω
voltage = 15V

# const outerRadius = 22.3mm |> m        #This governs how many layers of wires will be on the coil
clength = 1.0inch |> m
bthickness = 1.0mm
# const wireRadius = 1.6mm |> m      #This includes the insulation layer
numberOfCoils=4

coil = let
    innerRadius= projrad + bthickness
    outerRadius= innerRadius + 1.0inch |> m     #This governs how many layers of wires will be on the coil
    lengthCoil = 1.0inch |> m
    wireRadius = 1.6mm |> m       #This includes the insulation layer
    location = 0.0mm |> m 
    coilOnRange = 2.0*lengthCoil |> m 
    Coil(
    innerRadius,
    outerRadius,
    lengthCoil,
    wireRadius,
    location,
    coilOnRange
    )
end

barrel = let
    binnerRadius = projrad |> m
    thickness = bthickness |> m
    blength= numberOfCoils*clength |> m
    Barrel(
        binnerRadius,
        thickness,
        blength
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

#coils = CoilGenerator(numberOfCoils, projrad, projrad+barrel.thickness, length, wireRadius)


endTime = .2s

default_scenario = UIScenario(
    ip,
    barrel,
    coil,
    endTime,
    voltage,
    resistor,
    initialMagIRR,
    position,
    velocity,
    magnetization,
    numberOfCoils
)
function transform_scenario(scenario::UIScenario)::Scenario
    coils = CoilGenerator(scenario.numberOfCoils, scenario.coil.innerRadius, scenario.coil.outerRadius, scenario.coil.length, scenario.coil.wireRadius)
    PCE = ProjectileCoilEvent()
    PCE.entersActiveZone = [nothing for _ in coils]
    PCE.exitsActiveZone = [nothing for _ in coils]
    totalΩ = resistor + resistance(coils[1])
    scenario_to_be_solved = Scenario(
        scenario.ip,
        scenario.barrel,
        coils,
        PCE,
        scenario.endTime,
        scenario.voltage,
        scenario.resistor,
        scenario.initialMagIRR,
        scenario.initalPosition,
        scenario.initialVelocity,
        scenario.initialMagnetization,
    )
    return scenario_to_be_solved
end
function solve_scenario(scenario_to_be_solved::Scenario)
    solveScenario(scenario_to_be_solved)
end
end