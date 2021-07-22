module CoilGun
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn, 𝐈
using Unitful: Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Velocity, Time, Acceleration
import ForwardDiff
using Unitful

module CreatedUnits
using Unitful
using Unitful: 𝐈, 𝐌, 𝐓, 𝐋 , T, m, A, s
@derived_dimension HFieldGrad 𝐈*𝐋^-2
@derived_dimension HFieldRate 𝐈*𝐋^-1*𝐓^-1

@unit T_m "T/m" BFieldGradient 1T/m true
@unit A_ms "A/m/s" HFieldRate 1A/(m*s)      true
end


export IronProjectile
export NickelProjectile
export Coil
export Barrel
export ProjectilePhysical
export ProjectileMagnetic
export CoilGenerator
export ProjectileCoilEvent
export Scenario


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
    effectiveRange::Length   #range where the coil can affect the projectile
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
export volume
export mass
export density
export magDomainVol
export saturationMagnetizationFe
export numberWindings
export numberLayers
export totalNumberWindings
export wireLength
export area
export resistance
export coilCrossSectionalArea
export meanMagneticRadius
export acceleration
include("BasicFunctions.jl")


#Equations relating to the calculation of current
export ∂Current
export current
export ∂projectileInducedVoltage
export projectileInducedVoltage
export selfInductance
include("Current.jl")

# Functions for the magnetic field
export ∂HField_∂Current
export hFieldCoil
export ∇HFieldCoil
export dHField
include("MagneticField.jl")

# Functions for calculating material Magnetism
export δ
export δM
export ℒ
export ∂ℒ
export mag_Irr
export ∂Mag_irr_∂He
export ∂Magnetization_∂HField
include("Magnetization.jl")
#Force Functions

export totalForce
export dipoleCoilForce
export frictionForce
export airResistance
include("Forces.jl")

export solveScenario
export InitialConditions
include("solver.jl")

# include("MachineLearning.jl")

#export data to server
export dictionary_api
export get_default_scenario_json
include("api/api.jl")

end
#module

"""
Reminder: The pearmeability inside the coil is dependent upon the projectile's position, and the magnetic field from the coil acts differently inside the projectile

In order to maximize the force between the coil and the projectile, a moving sweet spot needs to be created. What this means is that the spot of highest gradient should stay at a constant distance from the projectile. The thinner the coils, the easier this can be accomplished. I would recommend investigating the optimal coil length in order to achieve this.
"""