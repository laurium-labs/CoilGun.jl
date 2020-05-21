module CoilGun

module CreatedUnits
    using Unitful
    using Unitful: 𝐈, 𝐌, 𝐓, 𝐋
    @derived_dimension BFieldGrad 𝐈^-1*𝐌*𝐓^-2*𝐋^-1
    @derived_dimension Permeability 𝐈/𝐋^2
    @derived_dimension HFieldRate 𝐈*𝐋^-1*𝐓^-1
end

using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn
using Unitful: Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Velocity, Time, Acceleration
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
const roomTemp = 300K                                       #Standard room Tempearture
const domainPinningFactor = 150A/m                          #This is the domain pinning factor for Iron (transformer).
const domainMagnetization = 0.2 * numberAtomsperDomainFe*bohrMagnetonPerAtomFe |> A/m #Magnetization of the domain
const magMomentPerDomain = domainMagnetization*domainSizeFe^3#This dipole magnetic moment doesn't take hysteresis/pinning into effect           #Unreasonably Small
const saturationMagnetizationPerKgFe = 217.6A/(m*kg)        #Saturation magnetizaiton of pure Iron per unit mass.
const kineticFrictionCoefficientFe = 0.36                   #Kinetic friction coefficient of Mild Steel on Copper, probably not exact
const staticFrictionCoefficientFe = 0.53                    #Static friction coefficient of copper on Steel, probably not exact
const dynamicViscosityAir = 1.825e-5kg/(m*s)                #Dynamic viscosity of air at 20C

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
    outerRadius::Length   #This governs how many layers of wires will be on the coil
    length::Length
    wireRadius::Length   #This includes the insulation layer
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
resistance(coil::Coil)::ElectricalResistance= resistivityCu*wireLength(coil)/area(coil)
coilCrossSectionalArea(coil::Coil)::Area    = (coil.outerRadius - coil.innerRadius) * coil.length
meanMagneticRadius(coil::Coil)::Length      = 2*coil.innerRadius*coil.outerRadius/(coil.innerRadius+coil.outerRadius)


#Equations relating to the calculation of current
function selfInductance(coil::Coil)::Inductance
    #This is a simplified version for the self inductance of the coil. It is not taking into consideration the thickness (different layers) of the coil, or the varying magnetic fields that pass through each loop. This is just for approximation only. When working out the math, current drops out of the equation so here it is just some random value.
    arbitraryCurrent = 1A
    return simpleBField(coil, arbitraryCurrent, 0m) * pi * totalNumberWindings(coil) * meanMagneticRadius(coil)^2/(3*arbitraryCurrent)
end
function projectileInducedVoltage(proj::Projectile, coil::Coil, magnetization::HField, velocity::Velocity, position::Length)::Voltage
    radius = meanMagneticRadius(coil)
    simpleArea = pi * radius^2
    ∂AreaRatio_∂t = radius * velocity * position/(position^2 + radius^2)^(3/2)
    constant = μ0 * magnetization * totalNumberWindings(coil) * simpleArea
    return constant * ∂AreaRatio_∂t
end
function ∂projectileInducedVoltage(coil::Coil, position::Length, velocity::Velocity, acceleration::Acceleration, magnetization::HField)
    #This funciton describes the change in the induced voltage per change in time
    radius = meanMagneticRadius(coil)
    simpleArea = pi * radius^2
    ∂∂AreaRatio_∂tt = radius * (position*acceleration/(position^2+radius^2)^(3/2) + radius^2*velocity^2/(position^2+radius^2)^(5/2)) |> s^-2
    constant = μ0 * magnetization * totalNumberWindings(coil) * simpleArea
    return constant * ∂∂AreaRatio_∂tt |> V/s
end
function coilCurrent(time, voltage, characteristicTime, couplingRelation, resistance)
    return voltage * (1-exp(-time / characteristicTime) * couplingRelation) / resistance
end
function current(proj::Projectile, coil::Coil, totalΩ::ElectricalResistance, voltage::Voltage, time::Time, magnetization::HField, velocity::Velocity, position::Length)::Current
    #This function calculates the current that is traveling through a coil. This is not taking operational amplifiers into consideration.
    arbitraryCurrent = 1A
    couplingFactor = simpleBField(coil, arbitraryCurrent, coil.length)/simpleBField(coil, arbitraryCurrent, 0m)
    constant = (1 - sqrt(1 - 4 * couplingFactor^2)) / (2 * couplingFactor^2)
    couplingRelation = exp(constant)/(1+constant*(constant-1))
    τ = selfInductance(coil)/totalΩ #Characteristic Time (When the current reaches 1-1/e of it's steady state value (5*τ))
    coilCurr(time::Time)::Current = coilCurrent(time,voltage,τ,couplingRelation,totalΩ)
    projectileInducedCurrent = projectileInducedVoltage(proj, coil, magnetization, velocity, position)/totalΩ
    return coilCurr(time) + projectileInducedCurrent
end
function ∂Current(coil::Coil, time::Time, voltage::Voltage, totalΩ::ElectricalResistance, position::Length, velocity::Velocity, acceleration::Acceleration, magnetization::HField)
    #This function describes how the current through the coil changes with the change in time.
    arbitraryCurrent = 1A
    couplingFactor = simpleBField(coil, arbitraryCurrent, coil.length)/simpleBField(coil, arbitraryCurrent, 0m)
    constant = (1 - sqrt(1 - 4 * couplingFactor^2)) / (2 * couplingFactor^2)
    couplingRelation = exp(constant)/(1+constant*(constant-1))
    τ = selfInductance(coil)/totalΩ |> s |> ustrip#Characteristic Time (When the current reaches 1-1/e of it's steady state value (5*τ))
    coilCurr(time) = coilCurrent(time,voltage |> V |> ustrip,τ,couplingRelation,totalΩ|> Ω |> ustrip)
    return ForwardDiff.derivative(coilCurr, time|>s |>ustrip)*1A/s + ∂projectileInducedVoltage(coil,position,velocity,acceleration,magnetization)/totalΩ
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
    return constant*mag(position)|> T
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
function ∂SimpleBField_∂Current(coil::Coil, current::Current, position::Length)
    #This function describes how the BField changes with respect to the change in current
    return simpleBField(coil, current, position)/current
end

#The paper referenced for these following equations relating to the magnetization of the projectile makes use of the Wiess mean Field theory in order to predict how the sample as a whole will react under a certain magnetic field.
δ(inc::CreatedUnits.HFieldRate)::Int = inc/sqrt(inc^2)
function δM(proj::Projectile, bField::BField, Mag_irr::HField, inc::CreatedUnits.HFieldRate)::Int
    #This corrects for when the field is reversed, and the difference between the irriversible magnetization (Mag_irr) and the and the anhysteris magnetization is the reversible magnetization. This function should take the values of 1 or 0.
    Mrev = proj.magnetic.saturationMagnetization * ℒ(proj, bField, Mag_irr) - Mag_irr
    dummyVar = Mrev/inc
    return (1 + dummyVar/sqrt(dummyVar^2))/2
end
function Mag_irr(proj::Projectile, bField::BField, Mag_irr::HField, magnetization::HField)::HField
    #This calculates the bulk irriversible magnetization inside the projectile.
    (magnetization - proj.magnetic.reversibility * ℒ(proj,bField,Mag_irr)*proj.magnetic.saturationMagnetization)/(1-proj.magnetic.reversibility)
end
function ∂Mag_irr_∂H(proj::Projectile, deltaM::Int, langevin::Float64, Mag_irr::HField)::HField
    return deltaM*(proj.magnetic.saturationMagnetization * langevin - Mag_irr)
end
function ℒ(proj::Projectile, bField::BField, Mag_irr::HField)::Float64
    #langevin funciton that represents the anhystesis bulk magnetization for a given material. It can be imagined as a sigmoid shape on a M-H graph.
    a = k*roomTemp/magMomentPerDomain |>T  |> ustrip               #Constant
    effectiveBField = bField+μ0*proj.magnetic.interdomainCoupling*Mag_irr |> T |> ustrip#Variable
    taylorApproxℒ = effectiveBField/(3*a) - effectiveBField^3/(45*a^3)
    return abs(effectiveBField/a) > 0.01 ? coth(effectiveBField/a) - a/effectiveBField : taylorApproxℒ
end
function ∂ℒ(proj::Projectile, bField::BField, Mag_irr::HField)::Float64
    #The first order derivative (with respect to the BField) of the ℒ function
    a = k*roomTemp/magMomentPerDomain |>T  |> ustrip               #Constant
    effectiveBField = bField+μ0*proj.magnetic.interdomainCoupling*Mag_irr |>T |>ustrip  #Variable
    ∂taylorApproxℒ = 1/(3*a) - effectiveBField^1/(15*a^3)
    langevin(x) = coth(x/a) - a/x
    return abs(effectiveBField/a) > 0.01 ? ForwardDiff.derivative(langevin,effectiveBField) : ∂taylorApproxℒ
end
function ∂HField(coil::Coil, current::Current, voltage::Voltage, totalΩ::ElectricalResistance,∇B::CreatedUnits.BFieldGrad, magnetization::HField, position::Length, velocity::Velocity, acceleration::Acceleration, time::Time)::CreatedUnits.HFieldRate
    #This function calculates the change in the HField due to the change in position and the change in current
    return (∇B*velocity+∂SimpleBField_∂Current(coil,current,position)*∂Current(coil,time,voltage,totalΩ,position,velocity,acceleration,magnetization))/μ0|>A/m/s
end
function ∂Magnetization_∂HField(proj::Projectile, bField::BField, Mag_irr::HField, ∂H::CreatedUnits.HFieldRate)::Float64
    #Change in the objects magnetization due to an external B-Field.
    ΔM_irr = (proj.magnetic.saturationMagnetization * ℒ(proj, bField,Mag_irr) - Mag_irr)
    numerator = δM(proj,bField,Mag_irr,∂H) * ΔM_irr + proj.magnetic.reversibility * ∂ℒ(proj, bField, Mag_irr) * domainPinningFactor*δ(∂H)
    denominator = (domainPinningFactor*δ(∂H)) - α * numerator
    return numerator/denominator
end
function ∂Magnetization(proj::Projectile, bField::BField, Mag_irr::HField, velocity::Velocity, ∇B::CreatedUnits.BFieldGrad, coil::Coil)::HField
    #Change in the objects magnetization due to an external B-Field.
    ΔH = (∇B * velocity + simpleBField(coil, I - prevI, position)/t) * t / μ0
    ΔM_irr = (proj.magnetic.saturationMagnetization * ℒ(proj, bField,Mag_irr) - Mag_irr)
    numerator = δM(proj,bField,Mag_irr,ΔH) * ΔM_irr + proj.magnetic.reversibility * ∂ℒ(proj, bField, Mag_irr) * (domainPinningFactor*δ(ΔH))
    denominator = (domainPinningFactor*δ(ΔH)) - α * numerator
    return ΔH * numerator/denominator
end

#Force Functions
function dipoleCoilForce(proj::Projectile, ∇BField::CreatedUnits.BFieldGrad, magnetization::HField)::Force
    #This force function assumes a number of things. 1) The magnetization of the projectile is constant throughout the material. This means that all of the magnetic domains are consistent throughout the material. 2) The magnetic field experienced at the center of the projectile is the average magnetic field experienced by the projectile. 3) The magnetization of the material can be approximated as a magnetic dipole (loop would be more accurate, but this is easier).
    magneticDipoleMoment = magnetization * volume(proj)
    return magneticDipoleMoment * ∇BField
end
function frictionForce(proj::Projectile, velocity::Velocity)::Force
    #This describes the resistive force due to friction (both static and kinetic). I know very little about the relationship between kinetic friction and velocity, but I highly doubt it's a constant relationship. More research will have to be done.
    frictionCoefficient = abs(velocity) > 0m/s ? kineticFrictionCoefficientFe : staticFrictionCoefficientFe
    normalForce = mass(proj) * gn
    return frictionCoefficient * normalForce
end
function airResistance(proj::Projectile, velocity::Velocity)::Force
    #This funciton is used to calculate the air resistance on the projectile. This current function is overly simplified and will need to be changed later for a more accurate function.
    return 6 * pi * dynamicViscosityAir * proj.physical.radius * -velocity
end
function totalForce(proj::Projectile, ∇BField::CreatedUnits.BFieldGrad, velocity::Velocity, magnetization::HField)::Force
    friction = frictionForce(proj, velocity)
    return abs(dipoleCoilForce(proj,∇BField, magnetization)) > abs(friction) ? dipoleCoilForce(proj,∇BField, magnetization) + friction + airResistance(proj, velocity) |> N : 0N
end

#Functions that calculate the change in velocity and change in position.
acceleration(force::Force, mass::Mass)::Acceleration = force/mass |>m/s^2
Δvel(acceleration::Acceleration, time::Time)::Velocity = acceleration * time |> m/s
Δpos(velocity::Velocity, time::Time)::Length = velocity * time |> m

include("solver.jl")
export IronProjectile, NickelProjectile, Coil, Barrel, volume, mass, density, numberWindings, numberLayers, 
    wireLength, area, volume, resistance, magDomainVol, magneticFieldSummation, magneticFieldIntegration, 
    ProjectilePhysical, ProjectileMagnetic, bFieldGradient,magDomainVol,saturationMagnetizationFe, coilCrossSectionalArea, 
    meanMagneticRadius, ℒ, ∂ℒ, dipoleCoilForce, totalNumberWindings, simpleBField, ∂Magnetization_∂HField, selfInductance, 
    projectileInducedVoltage, frictionForce, airResistance, current, totalForce, δ, δM , Mag_irr, ∂projectileInducedVoltage, 
    ∂Current, acceleration, ∂SimpleBField_∂Current, ∂HField, ∂Mag_irr_∂H,
    solveScenario, Scenario
end
#module

"""
Reminder: The pearmeability inside the coil is dependent upon the projectile's position, and the magnetic field from the coil acts differently inside the projectile

-Calc inductance of coil (use for time constant L/R)

Ferromagnetics doesn't have a linear susceptibility with an applied magnetic field and it's magnetization state if it isn't fully saturated. Meaning if the iron rod is fully saturated, the susceptibility is linear.

Ignoring:
domain rotation in k value for J-A model
"""