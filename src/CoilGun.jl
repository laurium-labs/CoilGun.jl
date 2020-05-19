module CoilGun

module CreatedUnits
    using Unitful
    using Unitful: 𝐈, 𝐌, 𝐓, 𝐋
    @derived_dimension BFieldGrad 𝐈^-1*𝐌*𝐓^-2*𝐋^-1
    @derived_dimension Permeability 𝐈/𝐋^2
end

using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn
using Unitful: Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Velocity, Time
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
mutable struct ProjectileMagnetic <: ElectroMagnetic
    domainSize::Length
    interdomainCoupling::Number
    magnetization::HField
    saturationMagnetization::HField
    reversibility::Float64
end
mutable struct IronProjectile <: Projectile
    physical::Physical
    magnetic::ElectroMagnetic
    position::Length       #This position is determined from the center of the coil to the center of the projectile
    velocity::Velocity
    #(I sense there'll be a problem when I try to incorperate multiple coils, but for now this is it).
end
mutable struct NickelProjectile <: Projectile
    physical::ProjectilePhysical
    magnetic::ProjectileMagnetic
    position::Length
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
function projectileInducedVoltage(proj::Projectile, coil::Coil)::Voltage
    radius = meanMagneticRadius(coil)
    simpleArea = pi * radius^2
    ∂AreaRatio_∂t = radius * proj.velocity * proj.position/(proj.position^2 + radius^2)^(3/2)
    constant = μ0 * proj.magnetic.magnetization * totalNumberWindings(coil) * simpleArea
    return constant * ∂AreaRatio_∂t
end
function current(proj::Projectile, coil::Coil, resistor::ElectricalResistance, voltage::Voltage, time::Time)::Current
    #This function calculates the current that is traveling through a coil. This is not taking operational amplifiers into consideration.
    totalΩ = resistor + resistance(coil)
    arbitraryCurrent = 1A
    couplingFactor = simpleBField(coil, arbitraryCurrent, coil.length)/simpleBField(coil, arbitraryCurrent, 0m)
    constant = (1 - sqrt(1 - 4 * couplingFactor^2)) / (2 * couplingFactor^2)
    couplingRelation = exp(constant)/(1+constant*(constant-1))
    τ = selfInductance(coil)/totalΩ #Characteristic Time (When the current reaches 1-1/e of it's steady state value (5*τ))
    coilCurrent = voltage * (1-exp(-time / τ) * couplingRelation) / totalΩ
    projectileInducedCurrent = projectileInducedVoltage(proj, coil)/totalΩ
    return coilCurrent + projectileInducedCurrent
end


# Funcitons for the magnetic field
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
    constant = μ0*totalNumberWindings(coil)*current/2
    mag(z::Length) = effectiveRadius^2/(effectiveRadius^2 + z^2)^(3/2)
    return constant*mag(position)
end
function simpleBField(coil::Coil, current::Current, coilPosition::Length, globalPosition::Length)::BField
    return simpleBField(coil, current, coilPosition-globalPosition)
end
function bFieldGradient(coil::Coil, current::Current, position::Length) :: CreatedUnits.BFieldGrad
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

#The paper referenced for these following equaitons relating to the magnetization of the projectile makes use of the Wiess mean Field theory in order to predict how the sample as a whole will react under a certain magnetic field.
δ(inc::HField)::Int = (inc > 0A/m) ? 1 : -1
function δM(proj::Projectile, bField::BField, Mag_irr::HField, inc::HField)::Int
    #This corrects for when the field is reversed, and the difference between the irriversible magnetization (Mag_irr) and the 
    return (proj.magnetic.saturationMagnetization * ℒ(proj, bField, Mag_irr) - Mag_irr)/inc >= 0 ? 1 : 0
end
function Mag_irr(proj::Projectile, bField::BField, Mag_irr::HField)::HField
    #This calculates the bulk irriversible magnetization inside the projectile.
    (proj.magnetic.magnetization - proj.magnetic.reversibility * ℒ(proj,bField,Mag_irr)*proj.magnetic.saturationMagnetization)/(1-proj.magnetic.reversibility)
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
function ΔMagnetization(proj::Projectile, bField::BField, Mag_irr::HField, ΔH::HField)::HField
    #Change in the objects magnetization due to an external B-Field.
    ΔM_irr = (proj.magnetic.saturationMagnetization * ℒ(proj, bField,Mag_irr) - Mag_irr)
    numerator = δM(proj,bField,Mag_irr,ΔH) * ΔM_irr + proj.magnetic.reversibility * ∂ℒ(proj, bField, Mag_irr) * (domainPinningFactor*δ(ΔH))
    denominator = (domainPinningFactor*δ(ΔH)) - α * numerator
    return ΔH * numerator/denominator
end

#Force Funcitons
function dipoleCoilForce(proj::Projectile, ∇BField::CreatedUnits.BFieldGrad)::Force
    #This force function assumes a number of things. 1) The magnetization of the projectile is constant throughout the material. This means that all of the magnetic domains are consistent throughout the material. 2) The magnetic field experienced at the center of the projectile is the average magnetic field experienced by the projectile. 3) The magnetization of the material can be approximated as a magnetic dipole (loop would be more accurate, but this is easier).
    magneticDipoleMoment = proj.magnetic.magnetization * volume(proj)
    return magneticDipoleMoment * ∇BField
end
function frictionForce(proj::Projectile)::Force
    #This describes the resistive force due to friction (both static and kinetic). I know very little about the relationship between kinetic friction and velocity, but I highly doubt it's a constant relationship. More research will have to be done.
    frictionCoefficient = abs(proj.velocity) > 0m/s ? kineticFrictionCoefficientFe : staticFrictionCoefficientFe
    normalForce = mass(proj) * gn
    return frictionCoefficient * normalForce
end
function airResistance(proj::Projectile)::Force
    #This funciton is used to calculate the air resistance on the projectile. This current function is overly simplified and will need to be changed later for a more accurate function.
    return 6 * pi * dynamicViscosityAir * proj.physical.radius * -proj.velocity
end
function totalForce(proj::Projectile, ∇BField::CreatedUnits.BFieldGrad)::Force
    return abs(dipoleCoilForce(proj,∇BField)) > abs(frictionForce(proj)) ? dipoleCoilForce(proj,∇BField) + frictionForce(proj) + airResistance(proj) |> N : 0N
end

#Functions that calculate the change in velocity and change in position.
function Δvel(proj::Projectile, force::Force, time::Time)::Velocity
    return force*time/mass(proj) |> m/s
end
function Δpos(proj::Projectile, time::Time)::Length
    return proj.velocity * time |> m
end


#Due to the nature of the changing magnetic field, creating an array of magnetic field is unreasonable because you would have to recalculate the array after each iteration.

# function generateBFieldGradient(coil::Coil,current::Current,proj::Projectile)
#     #This creates a gradient of the BField along the range of the BField.
#     function bFieldGradient(coil::Coil, current::Current, position::Length) :: CreatedUnits.BFieldGrad
#         effectiveRadius = meanMagneticRadius(coil) |> m |> ustrip
#         constant = μ0*totalNumberWindings(coil)*current/2
#         mag(z::Number)= effectiveRadius^2/(effectiveRadius^2 + z^2)^(3/2)
#         magGradient(z::Number) = ForwardDiff.derivative(mag, z)/m^2
#         position = position |> m |> ustrip
#         return constant*magGradient(position)
#     end 
#     integrationRange = 3*coil.length/2 |> m
#     dx = proj.magnetic.domainSize/2 |> m
#     coilPosition = 0m:dx:integrationRange
#     return [bFieldGradient(coil,current,x) for x in coilPosition]
# end
# function generateBField(coil::Coil,current::Current,proj::Projectile)
    # integrationRange = 3*coil.length/2 |> m
    # dx = proj.magnetic.domainSize/2 |> m
    # coilPosition = 0m:dx:integrationRange
    # return [simpleBField(x,coil,current) for x in coilPosition]
# end

# Magnetic Domain functions
# function generateMagneticDomians(physical::ProjectilePhysical, domainSize::Length, magneticStrengthperDomain::BField)

#     """A function that calculates the initial orientation of the magnetic domains along a 2-D slice of the projectile. """
#     numRings    = trunc(Int, physical.radius/domainSize)         #Number of concentric rings around the center of the rod that make up the domains (Rows)
#     numSlices   = trunc(Int, physical.length/domainSize)         #How many times the iron rod is sliced along the zAxis (Collumns)
#     #How the magnetic field and the domains interact
#     return [MagneticDipoleVector(exp(2*pi*rand()*im)*magneticStrengthperDomain) for slices in 1:numSlices, rings in 1:numRings]
# end
# function updateDomain(ironproj::Projectile,coil::Coil,bField::BField)
#     """
#     Here is where the magnetic moment of the domains is oriented with reference to the direction of the external magnetic field. I'm currently using trig to calculate the new direction using exponentials (e.g. exp[iθ]). This makes it easy to change the angle of the vector. This function seems pretty slow, and could probably be done more efficiently.
#     """
#     projLengthSize,radialLengthSize = size(ironproj.magnetic.domains)
#     coilEdgeToProjEdge = ironproj.position - (coil.length+ironproj.physical.length)/2
#     Δρ = ironproj.physical.radius/radialLengthSize
#     Δz =ironproj.physical.length/projLengthSize

#     #Both of these for loop probably could probably be replaced with a faster way to calculate the new domain orientation using mapping, but I know this works.
#     for ρ in 1:radialLengthSize #This for loop probably could be removed in place of creating an array with the desired values and calling them in the zAxis for loop
#         ρAxis = meanMagneticRadius(coil) - ρ*Δρ
#         for z in 1:projLengthSize
#             zAxis = z*Δz+coilEdgeToProjEdge
#             saturationAngle = atan(zAxis/ρAxis)+pi/2   #The MagneticDipoleVector must be pointed tangently from the meanMagneticRadius. There's probably a more efficient way to calculate this angle.
#             #What is being done here is i'm checking the orientation of the current magnetization, comparing it to the angle it would be at if the rod were fully saturated, and deciding which way from the saturation angle the new angle should be pointing.
#             δ = asin(imag(ironproj.magnetic.domains[z,ρ].vector))-saturationAngle > pi ? 1 : -1 #Ask Brent if calling the specific angle is correct
#             actualAngle = saturationAngle + acos((ironproj.magnetic.magnetization/saturationMagnetizationFe(ironproj)) |> ustrip) * δ
#             #This actual angle is accounting for domain wall movement and pinning and is based off of the B-H curve for magnetization of a ferromaterial.
#             θ = abs(zAxis) > coil.length/2 ? actualAngle : 0
#             ironproj.magnetic.domains[z,ρ] = MagneticDipoleVector(exp(θ*im)*ironproj.magnetic.magneticStrengthperDomain)
#         end
#     end
#     return ironproj.magnetic.domains
# end

# function magnetization(proj::Projectile, magField::BField, δ::Int)::HField
#     #This funciton is the basic funciton that the closing function and the effective magnetism is built out of. There are a couple different types: The normal funciton where theere are no special parameters, the reversal function that is the magnetization of the reversal point, the last magnetization which is the previous magnetization point, the + magnetization where the change in mag is positive, and correspondingly the - mag where the change is negative. This will have to be performed at each element of the projectile.
#     return proj.magnetic.saturationMagnetization * (langevin(proj, magField, 0)-domainPinningFactor*δ*langevin(proj, magField, 1)+domainPinningFactor^2*langevin(proj, magField, 2))
# end

# function closingFunction(magnetizationMinimum::HField, magnetizationMaximum::HField, proj::Projectile, bField::BField, prevBField::BField)::Float64
#     #Function insures that the B-H curve in the projectile has fixed endpoints and loops around that.
#     #The minimum and maximum magnetization points correspond to the max and min HFields that the point will experience.
#     if (bField - prevBField) >= 0T
#         δplus = 1
#         δminus = 0
#     else
#         δplus = 0
#         δminus = -1
#     end
#     return (magnetization(proj, bField, δminus)-magnetizationMinimum)/(magnetization(proj, bField, δplus)-magnetizationMaximum)
# end

# function effectiveMagnetization(proj::Projectile, bField::BField, bFieldMemory::Array{BField})::HField
#     #Question: Does the effective magnetizaiton describe the whole projectile, or could it describe a cell? No, using the Wiess Mean field theory this effective magnetization describes the entire projectile.
#     #magnetizationReverse is the point at where the HField reverses direction.
#     bMin = bFieldMemory[1]
#     bMax = bFieldMemory[2]
#     previousBField = bFieldMemory[3]
#     δ = (bField - previousBField) > 0T ? 1 : 0
#     magnetizationMinimum = magnetization(proj, bMin, -δ)
#     magnetizationMaximum = magnetization(proj, bMax, δ)
#     magnetizationReverse = (bField - previousBField) > 0T ? magnetizationMinimum : magnetizationMaximum
#     Λ = closingFunction(magnetizationMinimum, magnetizationMaximum, proj, bField, previousBField)
#     return Λ * (proj.magnetic.magnetization-magnetizationReverse)
# end
# #Push dealing with the whole projectile to the end, then sum the contributions of each cell. Meaning this is performed iteratively on each element.
# #Note: that the positonAlongProjectile is with respect to the coil. Only need to worry about force in X direction
# function domainCoilForce(positonAlongProjectile::Array{Int}, proj::Projectile, coil::Coil, bFieldGrad::CreatedUnits.BFieldGrad) :: Force
#     coilEdgeToProjEdge = proj.position - (coil.length+proj.physical.length)/2 |> m
#     heightDifference = meanMagneticRadius(coil) - positonAlongProjectile[2]*proj.magnetic.domainSize |> m
#     magneticMoment = proj.magnetic.magnetization * magDomainVol(proj) * real(proj.magnetic.domains[positonAlongProjectile[1],positonAlongProjectile[2]].angle)
#     return  magneticMoment * bFieldGrad * heightDifference/sqrt(heightDifference^2+coilEdgeToProjEdge^2)
# end

# function projectileCoilTotalForce(coil::Coil, proj::Projectile, ∇bField::BFieldGradient)
#     projLengthSize,radialLengthSize = size(proj.magnetic.domains)
#     elementSize = (3*coil.length/2)/size(∇bField.amplitude)[1] |> m #The span of the B-field

#     #The coordinateConversion function converts the incremental location within the projectile to the incremental location within the magnetic field.
#     coordinateConversion(x::Int) = Int(round((proj.position+(x-projLengthSize/2)*proj.magnetic.domainSize)/elementSize) |> ustrip)
#     totalForce = sum(dipoleCoilForce([z,ρ], proj, coil, ∇bField.amplitude[coordinateConversion(z),1])  for z = 1:projLengthSize for ρ = 1:radialLengthSize)
# end
#Is it benifitial to have matricies than arrays?
export IronProjectile, NickelProjectile, Coil, Barrel, volume, mass, density, numberWindings, numberLayers, wireLength, area, volume, resistance, magDomainVol, magneticFieldSummation, magneticFieldIntegration, MagneticDipoleVector, MagneticDipoleVector, ProjectilePhysical, ProjectileMagnetic, bFieldGradient,magDomainVol,saturationMagnetizationFe,coilCrossSectionalArea, meanMagneticRadius, generateBFieldGradient, generateMagneticDomians, updateDomain, ℒ, ∂ℒ, magnetization,closingFunction, effectiveMagnetization,dipoleCoilForce,projectileCoilTotalForce,totalNumberWindings, generateBField, simpleBField, ΔMagnetization, domainCoilForce, selfInductance, mutualInductance, projectileInducedVoltage, frictionForce, airResistance, current, totalForce, δ, δM , Mag_irr
end
#module

"""
Reminder: The pearmeability inside the coil is dependent upon the projectile's position, and the magnetic field from the coil acts differently inside the projectile

-Calc inductance of coil (use for time constant L/R)

Ferromagnetics doesn't have a linear susceptibility with an applied magnetic field and it's magnetization state if it isn't fully saturated. Meaning if the iron rod is fully saturated, the susceptibility is linear.

Ignoring:
domain rotation in k value for J-A model
"""