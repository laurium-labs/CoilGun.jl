module CoilGun

module CreatedUnits
  using Unitful
  using Unitful: 𝐈, 𝐌, 𝐓, 𝐋
  @derived_dimension BFieldGrad 𝐈^-1*𝐌*𝐓^-2*𝐋^-1
end
using Unitful:Ω, m, cm, kg, g, A, N, Na,ustrip, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm
using Unitful:Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance,BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density
using ForwardDiff
#MagneticDipoleMomentPerKg(::Unitful.FreeUnits{(A,m,kg), L^2*I*M^-1,nothing})

const resistivityCu = 1.72e-8m*Ω                            #Resistivity of Copper
const densityCu = 8960kg/m^3                                #Density if pure Copper
const densityFe = 7750kg/m^3                                #Density of pure Iron
const atomicWeightFe = 55.845g/mol                          #Atomic weight of Iron
const domainSizeFe = 26.5e-7m                               #Average magnetic domain size for pure Iron (actually 26.5e-9)
const densityNi = 8.908g/cm^3 |> kg/m^3                     #Density of pure Nickel
const currieTempFe = 1043K                                  #This is the Currie tempearture of Iron
const bohrMagnetonPerAtomFe = 2.2*μB                        #The is the Bohr Magneton per Iron atom
const numberAtomsperDomainFe = Na*domainSizeFe^3*densityFe/atomicWeightFe  #Number of atoms per Iron domain volume
const magPerFeAtom = currieTempFe*k/bohrMagnetonPerAtomFe   #This is the magnetic field given off by each Iron atom
const magPerFeDomain = magPerFeAtom*numberAtomsperDomainFe  #Magnetic field of the domain
const χFe = 200_000                                         #Magnetic susceptibility of iron at 20 C (unitless)
const μ = μ0*(1+χFe)                                        #Magnetic pearmeability of iron
const α = 9.5e-5                                            #Interdomain Coupling Factor (for an iron transformer)
const roomTemp = 300K                                       #Standard room Tempearture
const domainPinningFactor = 150                          #This is the domain pinning factor for Iron (transformer)
const domainMagnetization = numberAtomsperDomainFe*bohrMagnetonPerAtomFe  #Magnetization of the domain
const simplifiedMagMoment = domainMagnetization*domainSizeFe^3    #This dipole magnetic moment doesn't take hysteresis/pinning into effect
const saturationMagnetizationPerKgFe = 217.6A/(m*kg)             #Saturation magnetizaiton of pure Iron per unit mass.
#const meanMagneticRadius = 6.2438mm |> m                    #Radial position where the average magnetic field from the coils is located

abstract type Projectile end
abstract type Physical end
abstract type ElectroMagnetic end

struct MagneticDipoleVector
    angle::Complex{Float64}
    magnitude::BField #Using cylindrical coordindates
end

struct ProjectilePhysical <: Physical
    radius :: Length
    length :: Length
    density :: Density
end
mutable struct ProjectileMagnetic <: ElectroMagnetic
    domainSize::Length
    magneticStrengthperDomain::BField
    interdomainCoupling::Number
    magneticMoment::MagneticDipoleMoment
    magnetization::HField
    saturationMagnetization::HField
    domains :: Array{MagneticDipoleVector}
    magField :: BField
end
mutable struct IronProjectile <: Projectile
    physical::Physical
    magnetic::ElectroMagnetic
    position::Length       #This position is determined from the center of the coil to the center of the projectile
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
    thickness::Length   #This governs how many layers of wires will be on the coil
    length::Length
    wireRadius::Length   #This includes the insulation layer
end
struct BFieldGradient
    amplitude::Array{CreatedUnits.BFieldGrad}  #This should also include the position from the coil
end

#Below are functions associated with the projectile used
volume(proj::Projectile)::Volume    = proj.physical.radius^2 * π * proj.physical.length
mass(proj::Projectile)::Mass        = volume(proj)*density(proj)
density(proj::IronProjectile)       = proj.physical.density
magDomainVol(proj::Projectile)      = proj.magnetic.domainSize^3
saturationMagnetizationFe(proj::Projectile) = saturationMagnetizationPerKgFe*proj.physical.density*magDomainVol(proj)

#Below are functions associated with the Coil and wire
numberWindings(coil::Coil)          = trunc(Int,coil.length/(2*coil.wireRadius)) #Number of windings along the length of the coil
numberLayers(coil::Coil)            = trunc(Int,coil.thickness/(2*coil.wireRadius)) #Number of winding layers in the coil
totalNumberWindings(coil::Coil)     = trunc(Int,numberLayers(coil)*(numberWindings(coil)-0.5))
wireLength(coil::Coil)              = pi*numberLayers(coil)*(numberWindings(coil)*coil.innerRadius+sqrt(3)*coil.wireRadius*(numberWindings(coil)*(numberLayers(coil)+1)-(numberLayers(coil)+3)/2))
wireArea(coil::Coil) ::Area         = coil.wireRadius^2*pi
wireVolume(coil::Coil) ::Volume     = wireLength(coil)*wireArea(coil)
wireMass(coil::Coil) ::Mass         = densityCu*wireVolume(coil)
resistance(coil::Coil)              = resistivityCu*wireLength(coil)/wireArea(coil)
coilCrossSectionalArea(coil::Coil)  = (coil.innerRadius+coil.thickness) * coil.length
meanMagneticRadius(coil::Coil)      = 2*coil.innerRadius*(coil.thickness+coil.innerRadius)/(2*coil.innerRadius+coil.thickness)


#Below is the calculation of the magntic field from a coil through the summation of the BField from each individual loop.
function magneticFieldSummation(coil::Coil, current::Current, positionFromCoil::Length)
    coilRadius = coil.innerRadius
    wireRad = coil.wireRadius
    μ0/2*current*sum((layerNumber*sqrt(3)*wireRad+coilRadius)^2/((layerNumber*sqrt(3)*wireRad+coilRadius)^2+(positionFromCoil-(2*wireRad*(rowNumber-numberWindings(coil)/2)))^2)^(3/2) for rowNumber=1:numberWindings(coil) for layerNumber=1:numberLayers(coil)) 
end

#Reminder: The point starts in middle of the coil, then moves outward and goes through the front of the coil. When intPostion = coilLength it's at CoilFront.
function magneticFieldIntegration(coil::Coil, current::Current, integrationPositon::Length) :: BField
    coilInnerRadius = coil.innerRadius
    coilOuterRadius = coilInnerRadius + coil.thickness
    constant = μ0*totalNumberWindings(coil)*current/4
    function lengthIntegration(integrationPositon::Length)
        distToCoilBack = integrationPositon + coil.length/2
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
    return (constant * lengthIntegration(integrationPositon) / coilCrossSectionalArea(coil)) |> T
end
function magneticFieldIntegration(coil::Coil, current::Current, coilPosition::Length, globalPosition::Length) :: BField
    magneticFieldIntegration(coil, current, coilPosition - globalPosition)
end
function simpleBField(position::Length,coil::Coil, current::Current)::BField
    effectiveRadius = meanMagneticRadius(coil)
    constant = μ0*totalNumberWindings(coil)*current/2
    mag(z::Length) = effectiveRadius^2/(effectiveRadius^2 + z^2)^(3/2)
    return constant*mag(position)
end

#This needs to create a gradient of the BField along the integreation range.
function generateBFieldGradient(coil::Coil,current::Current,proj::Projectile)
    function bFieldGradient(coil::Coil, current::Current, position::Length)
        effectiveRadius = meanMagneticRadius(coil) |> m |> ustrip
        constant = μ0*totalNumberWindings(coil)*current/2
        mag(z::Number)= effectiveRadius^2/(effectiveRadius^2 + z^2)^(3/2)
        magGradient(z::Number) = ForwardDiff.derivative(mag, z)m
        position = position |> m |> ustrip
        return constant*magGradient(position)
    end 
    integrationRange = 3*coil.length/2 |> m
    dx = proj.magnetic.domainSize/2 |> m
    coilPosition = 0m:dx:integrationRange
    return [bFieldGradient(coil,current,x) for x in coilPosition]
end

function generateBField(coil::Coil,current::Current,proj::Projectile)
    integrationRange = 3*coil.length/2 |> m
    dx = proj.magnetic.domainSize/2 |> m
    coilPosition = 0m:dx:integrationRange
    return [simpleBField(x,coil,current) for x in coilPosition]
end

function generateMagneticDomians(ironproj::Projectile)
    """A function that calculates the initial orientation of the magnetic domains along a 2-D slice of the projectile. """
    numRings    = trunc(Int, ironproj.physical.radius/ironproj.magnetic.domainSize)         #Number of concentric rings around the center of the rod that make up the domains (Rows)
    numSlices   = trunc(Int, ironproj.physical.length/ironproj.magnetic.domainSize)         #How many times the iron rod is sliced along the zAxis (Collumns)
    #How the magnetic field and the domains interact
    ironproj.magnetic.domains = [MagneticDipoleVector(exp(2*pi*rand()*im),ironproj.magnetic.magneticStrengthperDomain) for slices in 1:numSlices, rings in 1:numRings]
 end

function updateDomain(ironproj::Projectile,coil::Coil,bField::BField)
    """
    Here is where the magnetic moment of the domains is oriented with reference to the direction of the external magnetic field. I'm currently using trig to calculate the new direction using exponentials (e.g. exp[iθ]). This makes it easy to change the angle of the vector. This function seems pretty slow, and could probably be done more efficiently.
    """
    projLengthSize,radialLengthSize = size(ironproj.magnetic.domains)
    coilEdgeToProjEdge = ironproj.position - (coil.length+ironproj.physical.length)/2
    Δρ = ironproj.physical.radius/radialLengthSize
    Δz =ironproj.physical.length/projLengthSize

    #Both of these for loop probably could probably be replaced with a faster way to calculate the new domain orientation using mapping, but I know this works.
    for ρ in 1:radialLengthSize #This for loop probably could be removed in place of creating an array with the desired values and calling them in the zAxis for loop
        ρAxis = meanMagneticRadius(coil) - ρ*Δρ
        for z in 1:projLengthSize
            zAxis = z*Δz+coilEdgeToProjEdge
            saturationAngle = atan(zAxis/ρAxis)+pi/2   #The MagneticDipoleVector must be pointed tangently from the meanMagneticRadius. There's probably a more efficient way to calculate this angle.
            #What is being done here is i'm checking the orientation of the current magnetization, comparing it to the angle it would be at if the rod were fully saturated, and deciding which way from the saturation angle the new angle should be pointing.
            δ = asin(imag(ironproj.magnetic.domains[z,ρ].angle))-saturationAngle > pi ? 1 : -1 #Ask Brent if calling the specific angle is correct
            actualAngle = saturationAngle + acos((ironproj.magnetic.magnetization/saturationMagnetizationFe(ironproj)) |> ustrip) * δ
            #This actual angle is accounting for domain wall movement and pinning and is based off of the B-H curve for magnetization of a ferromaterial.
            θ = abs(zAxis) > coil.length/2 ? actualAngle : 0
            ironproj.magnetic.domains[z,ρ] = MagneticDipoleVector(exp(θ*im),ironproj.magnetic.magneticStrengthperDomain)
        end
    end
    return ironproj.magnetic.domains
end

function langevin(proj::Projectile, bField::BField, derivative::Int64)::Float64
    a = k*roomTemp/proj.magnetic.magneticMoment |> T                #Constant
    effectiveBField = bField+μ0*α*proj.magnetic.magnetization |> T  #The effective B field is the B field experienced by the element
    x = effectiveBField/a |> ustrip                                     #Variable
    magnetization(var) = coth(var)-1/var
    if derivative > 0
        mag(var) = ForwardDiff.derivative(magnetization,var)
        if derivative > 1
            mag2(var) = ForwardDiff.derivative(var -> ForwardDiff.derivative(magnetization,var),var)
            return mag2(x)
        end
        return mag(x)
    end
    return magnetization(x)
end

function magnetization(proj::Projectile, magField::BField, δ::Int)
    #This funciton is the basic funciton that the closing function and the effective magnetism is built out of. There are a couple different types: The normal funciton where theere are no special parameters, the reversal function that is the magnetization of the reversal point, the last magnetization which is the previous magnetization point, the + magnetization where the change in mag is positive, and correspondingly the - mag where the change is negative. This will have to be performed at each element of the projectile.
    return proj.magnetic.saturationMagnetization * (langevin(proj, magField, 0)-domainPinningFactor*δ*langevin(proj, magField, 1)+domainPinningFactor^2*langevin(proj, magField, 2))
end
#########################################I'm really curious if the domainPinningFactor is unitless or not!

function closingFunction(magnetizationMinimum::HField, magnetizationMaximum::HField, proj::Projectile, bField::BField, prevBField::BField)::Float64
    #Function insures that the B-H curve in the projectile has fixed endpoints and loops around that.
    #The minimum and maximum magnetization points correspond to the max and min HFields that the point will experience.
    if (bField - prevBField) >= 0T
        δplus = 1
        δminus = 0
    else
        δplus = 0
        δminus = -1
    end
    return (magnetization(proj, bField, δminus)-magnetizationMinimum)/(magnetization(proj, bField, δplus)-magnetizationMaximum)
end

function effectiveMagnetization(proj::Projectile, bField::BField, bFieldMemory::Array{BField})::HField
    #Question: Does the effective magnetizaiton describe the whole projectile, or could it describe a cell?
    #magnetizationReverse is the point at where the HField reverses direction.
    bMin = bFieldMemory[1]
    bMax = bFieldMemory[2]
    magnetizationMinimum = langevin(proj,bMin,0)*saturationMagnetizationPerKgFe*mass(proj) #Fix
    magnetizationMaximum = langevin(proj,bMax,0)*saturationMagnetizationPerKgFe*mass(proj) #Fix
    previousBField = bFieldMemory[3]
    magnetizationReverse = (bField - previousBField) > 0T ? magnetizationMinimum : magnetizationMaximum
    Λ = closingFunction(magnetizationMinimum, magnetizationMaximum, proj, bField, previousBField)
    return Λ * (proj.magnetic.magnetization-magnetizationReverse)
end
#Push dealing with the whole projectile to the end, then sum the contributions of each cell. Meaning this is performed iteratively on each element.
#Note: that the positonAlongProjectile is with respect to the coil. Only need to worry about force in X direction
function dipoleCoilForce(positonAlongProjectile::Array{Int}, proj::Projectile, coil::Coil, bFieldGrad::CreatedUnits.BFieldGrad) :: Force
    coilEdgeToProjEdge = proj.position - (coil.length+proj.physical.length)/2 |> m
    heightDifference = meanMagneticRadius(coil) - positonAlongProjectile[2]*proj.magnetic.domainSize |> m
    magneticMoment = proj.magnetic.magnetization * magDomainVol(proj) * real(proj.magnetic.domains[positonAlongProjectile[1],positonAlongProjectile[2]].angle)
    return  magneticMoment * bFieldGrad * heightDifference/sqrt(heightDifference^2+coilEdgeToProjEdge^2)
end

function projectileCoilTotalForce(coil::Coil, proj::Projectile, ∇bField::BFieldGradient)
    projLengthSize,radialLengthSize = size(proj.magnetic.domains)
    elementSize = (3*coil.length/2)/size(∇bField.amplitude)[1] |> m #The span of the B-field

    #The coordinateConversion function converts the incremental location within the projectile to the incremental location within the magnetic field.
    coordinateConversion(x::Int) = Int(round((proj.position+(x-projLengthSize/2)*proj.magnetic.domainSize)/elementSize) |> ustrip)
    totalForce = sum(dipoleCoilForce([z,ρ], proj, coil, ∇bField.amplitude[coordinateConversion(z),1])  for z = 1:projLengthSize for ρ = 1:radialLengthSize)
end
#Is it benifitial to have matricies than arrays?
export IronProjectile, NickelProjectile, Coil, Barrel, volume, mass, density, numberWindings, numberLayers, wireLength, wireArea, wireVolume, wireMass, resistance, magDomainVol, magneticFieldSummation, magneticFieldIntegration, MagneticDipoleVector, MagneticDipoleVector, ProjectilePhysical, ProjectileMagnetic, BFieldGradient,magDomainVol,saturationMagnetizationFe,coilCrossSectionalArea, meanMagneticRadius, generateBFieldGradient, generateMagneticDomians, updateDomain,langevin,magnetization,closingFunction, effectiveMagnetization,dipoleCoilForce,projectileCoilTotalForce,totalNumberWindings, generateBField, simpleBField
end
#module

"""
Reminder: The pearmeability inside the coil is dependent upon the projectile's position, and the magnetic field from the coil acts differently inside the projectile

-Calc inductance of coil (use for time constant L/R)

Ferromagnetics doesn't have a linear susceptibility with an applied magnetic field and it's magnetization state if it isn't fully saturated. Meaning if the iron rod is fully saturated, the susceptibility is linear.

Ignoring:
domain rotation in k value for J-A model
"""