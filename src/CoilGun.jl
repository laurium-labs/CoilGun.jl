module CoilGun

module CreatedUnits
    using Unitful
    using Unitful: 𝐈, 𝐌, 𝐓, 𝐋 , T, m, A, s
    @derived_dimension BFieldGrad 𝐈^-1*𝐌*𝐓^-2*𝐋^-1
    @derived_dimension Permeability 𝐈/𝐋^2
    @derived_dimension HFieldRate 𝐈*𝐋^-1*𝐓^-1

    @unit T_m "T/m" BFieldGradient 1T/m true
    @unit A_ms "A/m/s" HFieldRate 1A/(m*s)      true
end

using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn, 𝐈
using Unitful: Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Velocity, Time, Acceleration
using ForwardDiff
using Unitful
using DifferentialEquations


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

#Magnetism Equation Parameters
const domainPinningFactor = 742.64A/m           #This is the domain pinning factor from Ref.[5]
const α = 1.34e-3                               #Interdomain Coupling Factor from Ref.[5]
const a = 882.55A/m                             #"Determines the density distribution of mag. domians"~Ref.[2] Ref.[5]
const magMomentPerDomain = k*roomTemp/(a * μ0)  #This dipole magnetic moment from Ref.[5]

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
function Coil(numberOfCoils::Int, innerRadius::Length, outerRadius::Length, coilLength::Length, wireRadius::Length)
    return[Coil(innerRadius, outerRadius, coilLength, wireRadius, x*coilLength , 2*coilLength) for x in 1:numberOfCoils]
end

mutable struct ProjectileCoilEvent
    entersActiveZone::Array{Union{Nothing,Time},1}
    exitsActiveZone::Array{Union{Nothing,Time},1}
    function ProjectileCoilEvent()
        new([nothing],[nothing])
    end
end

#Below are functions associated with the projectile used
volume(proj::Projectile)::Volume            = proj.physical.radius^2 * π * proj.physical.length
mass(proj::Projectile)::Mass                = volume(proj)*density(proj)
density(proj::Projectile)::Density          = proj.physical.density
magDomainVol(proj::Projectile)::Volume      = proj.magnetic.domainSize^3
saturationMagnetizationFe(proj::Projectile)::HField = saturationMagnetizationPerKgFe*proj.physical.density*volume(proj)

#Below are functions associated with the Coil and wire
numberWindings(coil::Coil)::Int             = trunc(Int,coil.length/(2*coil.wireRadius)) #Number of windings along the length of the coil
numberLayers(coil::Coil)::Int               = trunc(Int,(coil.outerRadius-coil.innerRadius)/(2*coil.wireRadius)) #Number of winding layers in the coil
totalNumberWindings(coil::Coil)::Int        = trunc(Int,numberLayers(coil)*(numberWindings(coil)-0.5))
wireLength(coil::Coil)::Length              = pi*numberLayers(coil)*(numberWindings(coil)*coil.innerRadius+sqrt(3)*coil.wireRadius*(numberWindings(coil)*(numberLayers(coil)+1)-(numberLayers(coil)+3)/2))
area(coil::Coil)::Area                      = coil.wireRadius^2*pi
volume(coil::Coil)::Volume                  = wireLength(coil)*area(coil)
mass(coil::Coil)::Mass                      = densityCu*volume(coil)
resistance(coil::Coil)::ElectricalResistance= resistivityCu*wireLength(coil)/area(coil) |> Ω
coilCrossSectionalArea(coil::Coil)::Area    = (coil.outerRadius - coil.innerRadius) * coil.length
meanMagneticRadius(coil::Coil)::Length      = 2*coil.innerRadius*coil.outerRadius/(coil.innerRadius+coil.outerRadius)


#Equations relating to the calculation of current
function selfInductance(coil::Coil)::Inductance
    #This is a simplified version for the self inductance of the coil. It is not taking into consideration the thickness (different layers) of the coil, or the varying magnetic fields that pass through each loop. This is just for approximation only. When working out the math, current drops out of the equation so here it is just some random value.
    iR = coil.innerRadius
    oR = coil.outerRadius
    length = coil.length
    B_I = 4*pi*μ0*totalNumberWindings(coil)/(9*length^2*(oR-iR)^2)
    var = ((oR^2+length^2)^(3/2)-(iR^2+length^2)^(3/2)-(oR^3-iR^3))*(oR^2-iR^2)|> m^5
    return B_I*var
end
function projectileInducedVoltage(coil::Coil, magnetization::HField, velocity::Velocity, position::Length)::Voltage #Fix
    radius = meanMagneticRadius(coil)
    position = position - coil.location
    simpleArea = pi * radius^2
    ∂AreaRatio_∂t = radius * velocity * position/(position^2 + radius^2)^(3/2)
    constant = μ0 * magnetization * totalNumberWindings(coil) * simpleArea
    return constant * ∂AreaRatio_∂t
end
function ∂projectileInducedVoltage(coil::Coil, position::Length, velocity::Velocity, acceleration::Acceleration, magnetization::HField)
    #This funciton describes the change in the induced voltage per change in time
    position = position - coil.location
    radius = meanMagneticRadius(coil)
    simpleArea = pi * radius^2
    ∂∂AreaRatio_∂tt = radius * (position*acceleration/(position^2+radius^2)^(3/2) + radius^2*velocity^2/(position^2+radius^2)^(5/2)) |> s^-2
    constant = μ0 * magnetization * totalNumberWindings(coil) * simpleArea
    return constant * ∂∂AreaRatio_∂tt |> V/s
end
function couplingFactor(coil::Coil)::Float64 #May need to be changed. There's another equation with number of turns
    #This function calculates the ratio of magnetic field lines passing through the generating coil, and an adjacent coil.
    a = coil.length/2
    α1 = 4 * a
    β1 = 2 * a
    β2 = 0m
    block(var::Length) = (coil.outerRadius^2 + var^2)^(3/2) - (coil.innerRadius^2 + var^2)^(3/2)
    numerator = block(α1) - 2*block(β1) + block(β2)
    denominator = 2 * (block(β1)-block(β2))
    return numerator/denominator
end

#TODO: Create a function that describes how a capacitor will supply voltage to a coil
#Functions for current
function coilCurrent(coil::Coil, position::Length, time::Time, maxVoltage::Voltage, characteristicTime::Time, resistance::ElectricalResistance)::Current #Fix: Turn coil off when proj. is > 2 coilOnRange away from coil
    #The entry fields are intentionally left without specification due to its derivative being taken.
    distFromCoil = coil.location - position
    if (0m < distFromCoil) && (distFromCoil <= coil.coilOnRange)
        return maxVoltage * (1-exp(-time / characteristicTime)) / resistance
    elseif distFromCoil <= 0m
        return maxVoltage * (2*exp(-time / characteristicTime)-1) / resistance
    else
        return maxVoltage/resistance - maxVoltage/resistance
    end
end
function current(coil::Coil, totalΩ::ElectricalResistance, initialVoltage::Voltage, time::Time, magnetization::HField, velocity::Velocity, position::Length)::Current
    #This function calculates the current that is traveling through a coil. This is not taking operational amplifiers into consideration.
    𝓀 = couplingFactor(coil)
    τ = selfInductance(coil)*𝓀^2/totalΩ #Characteristic Time (current reaches steady state value at 5*τ)
    projectileInducedCurrent = projectileInducedVoltage(coil, magnetization, velocity, position)/totalΩ
    return coilCurrent(coil, position|>m, time,initialVoltage,τ,totalΩ)|>A #+ projectileInducedCurrent |> A
end
function ∂Current(coil::Coil, time::Time, initialVoltage::Voltage, totalΩ::ElectricalResistance, position::Length)
    #This function describes how the current through the coil changes with the change in time.
    𝓀 = couplingFactor(coil)
    τ = selfInductance(coil)*𝓀^2/totalΩ #Characteristic Time (Current reaches it's steady state value at 5*τ)
    distFromCoil = coil.location - position
    ∂CoilVoltageRate = if (0m < distFromCoil) && (distFromCoil <= coil.coilOnRange)
        initialVoltage * exp(-time / τ)/τ
    elseif distFromCoil <= 0m
        -2*initialVoltage * exp(-time / τ)/τ
    else
        0V/s
    end
    return ∂CoilVoltageRate/totalΩ|>A/s# + ∂projectileInducedVoltage(coil,position,velocity,acceleration,magnetization))/totalΩ |> A/s
end

# Functions for the magnetic field
#Reminder: The point starts in middle of the coil, then moves outward and goes through the front of the coil. When intPostion = coilLength it's at CoilFront.
function magneticFieldSummation(coil::Coil, current::Current, positionFromCoil::Length)::BField
    #This calculates the magntic field from a coil through the summation of the BField from each individual loop.
    coilRadius = coil.innerRadius
    wireRad = coil.wireRadius
    μ0/2*current*sum((layerNumber*sqrt(3)*wireRad+coilRadius)^2/((layerNumber*sqrt(3)*wireRad+coilRadius)^2+(positionFromCoil-(2*wireRad*(rowNumber-numberWindings(coil)/2)))^2)^(3/2) for rowNumber=1:numberWindings(coil) for layerNumber=1:numberLayers(coil)) 
end 
function magneticFieldIntegration(coil::Coil, current::Current, positon::Length) :: BField
    coilInnerRadius = coil.innerRadius
    coilOuterRadius = coil.outerRadius
    constant = μ0*totalNumberWindings(coil)*current/4
    function lengthIntegration(positon::Length)
        distToCoilBack = positon + coil.length/2
        distToCoilFront = distToCoilBack - coil.length
        logarithm(coilPosition::Length,radialLength::Length) = log((coilPosition |> m |> ustrip)^2+(radialLength |> m |> ustrip)^2)
        arctan(coilPosition::Length,radialLength::Length) = atan(((coilPosition |> m) / (radialLength |> m))|>ustrip)
        ∫_radial(radialLength::Length) = 
            distToCoilBack  * logarithm(distToCoilBack,radialLength) - 
            distToCoilFront * logarithm(distToCoilFront,radialLength) +
            (distToCoilBack-distToCoilFront) + 2 * radialLength * 
            (arctan(distToCoilBack,radialLength) - 
            arctan(distToCoilFront,radialLength))
        return ∫_radial(coilOuterRadius) - ∫_radial(coilInnerRadius)
    end
    return (constant * lengthIntegration(positon) / coilCrossSectionalArea(coil)) |> T
end
function magneticFieldIntegration(coil::Coil, current::Current, coilPosition::Length, globalPosition::Length) :: BField
    magneticFieldIntegration(coil, current, coilPosition - globalPosition)
end
function simpleBField(coil::Coil, current::Current, position::Length)::BField
    effectiveRadius = meanMagneticRadius(coil)
    constant = μ0*totalNumberWindings(coil)*(current)/2
    mag(z) = effectiveRadius^2/(effectiveRadius^2 + z^2)^(3/2)
    return constant*mag(coil.location - position)|> T
end
function simpleBField(coil::Coil, current::Current, coilPosition::Length, globalPosition::Length)::BField
    return simpleBField(coil, current, coilPosition-globalPosition)
end
function bFieldGradient(coil::Coil, current::Current, position::Length)::CreatedUnits.BFieldGrad
    effectiveRadius = meanMagneticRadius(coil) |> m |> ustrip
    constant = μ0*totalNumberWindings(coil)*current/2
    mag(z::Number)= effectiveRadius^2/(effectiveRadius^2 + z^2)^(3/2)
    magGradient(z::Number) = ForwardDiff.derivative(mag, z)/m^2
    position = position |> m |> ustrip
    return constant*magGradient(position)
end 
function bFieldGradient(coil::Coil, current::Current, coilPosition::Length, globalPosition::Length)::CreatedUnits.BFieldGrad
    return bFieldGradient(coil, current, coilPosition-globalPosition)
end
function ∂BField_∂Current(coil::Coil, position::Length)
    #This function describes how the BField changes with respect to the change in current
    return bFieldCoil(coil, 1A, position)/1A
end
# function ∇BFieldCoil(coil::Coil, current::Current, position::Length)::CreatedUnits.BFieldGrad
#     # This function describs the change in the magnetic field of a coil with the change in position
#     innerRad = coil.innerRadius
#     outerRad = coil.outerRadius
#     length = coil.length
#     crossSectionalArea = (outerRad - innerRad)*length
#     constant = (μ0/(crossSectionalArea*(outerRad - innerRad)))
#     farEdgeofCoil = coil.location - position + length/2
#     closeEdgeofCoil = farEdgeofCoil - length
#     variable(a) = a*(sqrt(a^2+(outerRad|>ustrip)^2)-sqrt(a^2 + (innerRad|>ustrip)^2))
#     ∇variable(a::Length)::Length = -ForwardDiff.derivative(variable, a|>ustrip)m
#     return constant*current*totalNumberWindings(coil)*(∇variable(farEdgeofCoil)-∇variable(closeEdgeofCoil)) |> T/m
# end
function bFieldCoil(coil::Coil, current::Current, position::Length)::BField
    constant = μ0*current*totalNumberWindings(coil)/coilCrossSectionalArea(coil)
    logarithm(pos::Length)::Length = pos * log((sqrt(pos^2+coil.outerRadius^2)+coil.outerRadius)/(sqrt(pos^2+coil.innerRadius^2)+coil.innerRadius))
    farEdgeofCoil = coil.location - position + coil.length/2
    closeEdgeofCoil = farEdgeofCoil - coil.length
    return constant * (logarithm(farEdgeofCoil) - logarithm(closeEdgeofCoil))
end
function ∇BFieldCoil(coil::Coil, current::Current, position::Length)::CreatedUnits.BFieldGrad
    constant = μ0*current*totalNumberWindings(coil)/coilCrossSectionalArea(coil)
    outerRadius = coil.outerRadius |> m |>ustrip
    innerRadius = coil.innerRadius |>m |> ustrip
    logarithm(pos) = pos * log((sqrt(pos^2+outerRadius^2)+outerRadius)/(sqrt(pos^2+innerRadius^2)+innerRadius))
    ∇logarithm(pos::Length)::Float64 = -ForwardDiff.derivative(logarithm,pos|>m|>ustrip)
    farEdgeofCoil = coil.location - position + coil.length/2
    closeEdgeofCoil = farEdgeofCoil - coil.length
    return constant * (∇logarithm(farEdgeofCoil) - ∇logarithm(closeEdgeofCoil))
end

#The paper referenced for these following equations relating to the magnetization of the projectile makes use of the Wiess mean Field theory in order to predict how the sample as a whole will react under a certain magnetic field.
function δ(inc::CreatedUnits.HFieldRate)::Int
    return inc >= 0A/m/s ? 1 : -1
end
function δM(proj::Projectile, bField::BField, Mag_irr::HField, inc::CreatedUnits.HFieldRate)::Int
    #This corrects for when the field is reversed, and the difference between the irriversible magnetization (Mag_irr) and the and the anhysteris magnetization is the reversible magnetization. This function should take the values of 1 or 0.
    Mrev = ℒ(proj, bField, Mag_irr) - Mag_irr
    dummyVar = abs(Mrev) < 1e-6A/m ? 1 : Mrev/δ(inc)
    return (1 + sign(dummyVar))/2
end
function ℒ(proj::Projectile, bField::BField, Mag_irr::HField)::HField
    #langevin funciton that represents the anhystesis bulk magnetization for a given material. It can be imagined as a sigmoid shape on a M-H graph.
    a = k*roomTemp/magMomentPerDomain |>T  |> ustrip               #Constant
    effectiveBField = bField+μ0*proj.magnetic.interdomainCoupling*Mag_irr |> T |> ustrip#Variable
    taylorApproxℒ = effectiveBField/(3*a) - effectiveBField^3/(45*a^3)
    ans = abs(effectiveBField/a) > 0.01 ? coth(effectiveBField/a) - a/effectiveBField : taylorApproxℒ
    return ans * proj.magnetic.saturationMagnetization
end
function ∂ℒ(proj::Projectile, bField::BField, Mag_irr::HField)::Float64
    #The first order derivative (with respect to the BField) of the ℒ function
    a = k*roomTemp/magMomentPerDomain |>T|>ustrip             #Constant
    effectiveBField = bField+μ0*proj.magnetic.interdomainCoupling*Mag_irr |>T|>ustrip  #Variable
    ∂taylorApproxℒ = 1/(3*a) - effectiveBField^2/(15*a^3)
    langevin(x) = coth(x/a) - a/x
    ans = abs(effectiveBField/a) > 1e-6 ? ForwardDiff.derivative(langevin,effectiveBField)/1T : ∂taylorApproxℒ/1T
    # println("∂ℒ :\t Mag_irr $(μ0*proj.magnetic.interdomainCoupling*Mag_irr),\tRatio: $(effectiveBField),\tans: $(ans),\ta: $(a)")
    return ans * proj.magnetic.saturationMagnetization * μ0
end
function mag_Irr(proj::Projectile, bField::BField, Mag_irr::HField, magnetization::HField)::HField
    #This calculates the bulk irriversible magnetization inside the projectile.
    return (magnetization - proj.magnetic.reversibility * ℒ(proj,bField,Mag_irr))/(1-proj.magnetic.reversibility)
end
function ∂Mag_irr_∂He(proj::Projectile, bField::BField, Mag_irr::HField, ∂H::CreatedUnits.HFieldRate)::Float64
    return δM(proj,bField,Mag_irr,∂H)*(ℒ(proj,bField,Mag_irr) - Mag_irr)/(domainPinningFactor)
end
function dHField(coils::Array{Coil,1}, voltage::Voltage, totalΩ::ElectricalResistance, ∇B::CreatedUnits.BFieldGrad, position::Length, velocity::Velocity, time::Time)::CreatedUnits.HFieldRate
    #This function calculates the change in the HField due to the change in position and the change in current
    return (∇B*velocity+sum(map(coil -> ∂BField_∂Current(coil,position)*∂Current(coil,time,voltage,totalΩ,position), coils)))/μ0|>A/m/s
end

#Somehow the rod is oversaturating
function ∂Magnetization_∂HField(proj::Projectile, bField::BField, Mag_irr::HField, ∂H::CreatedUnits.HFieldRate)::Float64
    #Change in the objects magnetization due to an external B-Field.
    ΔM_irr = ∂Mag_irr_∂He(proj, bField, Mag_irr, ∂H)
    numerator = ΔM_irr + proj.magnetic.reversibility * ∂ℒ(proj, bField, Mag_irr) * δ(∂H)
    # println("∂Magnetization_∂HField:\tδ: $(δ(∂H))")
    denominator = δ(∂H) - α * numerator
    return numerator/denominator
end

#Force Functions
function dipoleCoilForce(proj::Projectile, ∇BField::CreatedUnits.BFieldGrad, magnetization::HField)::Force #Update: Integrate mag over volume
    #This force function assumes a number of things. 1) The magnetization of the projectile is constant throughout the material. This means that all of the magnetic domains are consistent throughout the material. 2) The magnetic field experienced at the center of the projectile is the average magnetic field experienced by the projectile. 3) The magnetization of the material can be approximated as a magnetic dipole (loop would be more accurate, but this is easier).
    magneticDipoleMoment = magnetization * volume(proj)
    return magneticDipoleMoment * ∇BField|>N
end
function frictionForce(proj::Projectile, velocity::Velocity, dipoleCoilForce::Force)::Force #Update (check inside)
    #This describes the resistive force due to friction (both static and kinetic). I know very little about the relationship between kinetic friction and velocity, but I highly doubt it's a constant relationship. More research will have to be done.
    frictionCoefficient = abs(velocity) > 1e-5m/s ? kineticFrictionCoefficientFe : staticFrictionCoefficientFe
    normalForce = mass(proj) * gn
    signCorrection = velocity == 0 ? 1 : -sign(velocity)
    if abs(normalForce) >= abs(dipoleCoilForce) && abs(velocity) < 1e-5m/s
        return -dipoleCoilForce|>N
    else
        return frictionCoefficient * normalForce * signCorrection|>N
    end
end
function airResistance(proj::Projectile, velocity::Velocity)::Force #Update to model rod through tube
    #This funciton is used to calculate the air resistance on the projectile. This current function is overly simplified and will need to be changed later for a more accurate function.
    return 6 * pi * dynamicViscosityAir * proj.physical.radius * velocity|>N
end
function totalForce(proj::Projectile, ∇BField::CreatedUnits.BFieldGrad, velocity::Velocity, magnetization::HField)::Force
    dipoleForce = dipoleCoilForce(proj,∇BField, magnetization)
    frictionalForces = frictionForce(proj, velocity, dipoleForce) + airResistance(proj, velocity)
    return dipoleForce + frictionalForces |> N
end
acceleration(force::Force, mass::Mass)::Acceleration = force/mass |>m/s^2

include("solver.jl")
export IronProjectile, NickelProjectile, Coil, Barrel, volume, mass, density, numberWindings, numberLayers, 
    wireLength, area, volume, resistance, magDomainVol, magneticFieldSummation, magneticFieldIntegration, 
    ProjectilePhysical, ProjectileMagnetic, bFieldGradient,magDomainVol,saturationMagnetizationFe, coilCrossSectionalArea, 
    meanMagneticRadius, ℒ, ∂ℒ, dipoleCoilForce, totalNumberWindings, simpleBField, ∂Magnetization_∂HField, selfInductance, 
    projectileInducedVoltage, frictionForce, airResistance, current, totalForce, δ, δM , Mag_irr, ∂projectileInducedVoltage, 
    ∂Current, acceleration, ∂BField_∂Current, dHField, bFieldCoil, ∇BFieldCoil, ProjectileCoilEvent,
    solveScenario, Scenario
end
#module

"""
Reminder: The pearmeability inside the coil is dependent upon the projectile's position, and the magnetic field from the coil acts differently inside the projectile

In order to maximize the force between the coil and the projectile, a moving sweet spot needs to be created. What this means is that the spot of highest gradient should stay at a constant distance from the projectile. The thinner the coils, the easier this can be accomplished. I would recommend investigating the optimal coil length in order to achieve this.
"""