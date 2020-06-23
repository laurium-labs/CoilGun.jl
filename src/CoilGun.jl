module CoilGun

module CreatedUnits
    using Unitful
    using Unitful: 𝐈, 𝐌, 𝐓, 𝐋 , T, m, A, s
    @derived_dimension HFieldGrad 𝐈*𝐋^-2
    @derived_dimension HFieldRate 𝐈*𝐋^-1*𝐓^-1

    @unit T_m "T/m" BFieldGradient 1T/m true
    @unit A_ms "A/m/s" HFieldRate 1A/(m*s)      true
end

using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn, 𝐈
using Unitful: Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Velocity, Time, Acceleration
using ForwardDiff
using Unitful
using DifferentialEquations


include("Constants.jl")
#Magnetism Equation Parameters
const domainPinningFactor = 742.64A/m           #This is the domain pinning factor from Ref.[5]
const α = 1.34e-3                               #Interdomain Coupling Factor from Ref.[5]
const a = 882.55A/m                             #"Determines the density distribution of mag. domians"~Ref.[2] Ref.[5]
const magMomentPerDomain = k*roomTemp/a  #This dipole magnetic moment from Ref.[5]

abstract type Projectile end
abstract type Physical end
abstract type ElectroMagnetic end

struct ProjectilePhysical <: Physical
    radius :: Length
    length :: Length
    density :: Density
end
struct ProjectileMagnetic <: ElectroMagnetic
    domainSize::Length
    interdomainCoupling::Number
    saturationMagnetization::HField
    reversibility::Float64
end
struct IronProjectile <: Projectile
    physical::Physical
    magnetic::ElectroMagnetic
end
struct NickelProjectile <: Projectile
    physical::ProjectilePhysical
    magnetic::ProjectileMagnetic
end
struct Barrel
    innerRadius::Length
    thickness::Length
    length::Length
end
struct Coil
    innerRadius::Length
    outerRadius::Length     #This governs how many layers of wires will be on the coil
    length::Length
    wireRadius::Length      #This includes the insulation layer
    location::Length        #Global position of the coil
    coilOnRange::Length   #section of barrel were coil is on
end
function CoilGenerator(numberOfCoils::Int, innerRadius::Length, outerRadius::Length, coilLength::Length, wireRadius::Length)
    return[Coil(innerRadius, outerRadius, coilLength, wireRadius, x*coilLength , 2*coilLength) for x in 1:numberOfCoils]
end

mutable struct ProjectileCoilEvent
    entersActiveZone::Array{Union{Nothing,Time},1}
    exitsActiveZone::Array{Union{Nothing,Time},1}
    function ProjectileCoilEvent()
        new([nothing],[nothing])
    end
end

#Functions relating to simple calculations:
include("BasicFunctions.jl")

include("api/api.jl")

#Equations relating to the calculation of current
include("Current.jl")
# Functions for the magnetic field
include("MagneticField.jl")
# Functions for calculating material Magnetism
include("Magnetization.jl")
#Force Functions
include("Forces.jl")
acceleration(force::Force, mass::Mass)::Acceleration = force/mass |>m/s^2

include("solver.jl")
export IronProjectile, NickelProjectile, Coil, Barrel, volume, mass, density, numberWindings, numberLayers, 
    wireLength, area, volume, resistance, magDomainVol, ∂Mag_irr_∂He,ProjectilePhysical, ProjectileMagnetic, 
    magDomainVol,saturationMagnetizationFe, coilCrossSectionalArea, CoilGenerator,
    meanMagneticRadius, ℒ, ∂ℒ, dipoleCoilForce, totalNumberWindings, ∂Magnetization_∂HField, selfInductance, 
    projectileInducedVoltage, frictionForce, airResistance, current, totalForce, δ, δM , mag_Irr, ∂projectileInducedVoltage, 
    ∂Current, acceleration, ∂HField_∂Current, dHField, hFieldCoil, ∇HFieldCoil, ProjectileCoilEvent,
    solveScenario, Scenario, export_array
end
#module

"""
Reminder: The pearmeability inside the coil is dependent upon the projectile's position, and the magnetic field from the coil acts differently inside the projectile

In order to maximize the force between the coil and the projectile, a moving sweet spot needs to be created. What this means is that the spot ofç highest gradient should stay at a constant distance from the projectile. The thinner the coils, the easier this can be accomplished. I would recommend investigating the optimal coil length in order to achieve this.
"""