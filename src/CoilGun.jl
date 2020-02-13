module CoilGun
using Unitful:Ω, m, cm, kg, A, N, ustrip, T, s, μ0, ϵ0, Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area

const resistivityCu = 1.72e-8m*Ω    #Resistivity of Copper
const densityCu = 8960kg/m^3        #Density if pure Copper
const densityFe = 7750kg/m^3        #Density of pure Iron
const domainSizeFe = 26.5e-9m    #Average magnetic domain size for pure Iron
const densityNi = 7750kg/m^3 # fix
abstract type Projectile end
struct IronProjectile <: Projectile
 radius::Length
 length::Length
 domainSize::Length
end
struct NickelProjectile <: Projectile
    radius::Length
    length::Length
    domainSize::Length
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
    wireRadius::Length   #This is the wires radius including the insulation layer
end

#Below are functions associated with the projectile used
volume(proj::Projectile)::Volume = proj.radius^2 * π * proj.length
mass(proj::Projectile)::Mass = volume(proj)*density(proj)
density(proj::IronProjectile) = densityFe
density(proj::NickelProjectile) = densityNi 

#Below are functions associated with the Coil and wire
numberWindings(coil::Coil) = trunc(Int,coil.length/(2*coil.wireRadius)) #Number of windings along the length of the coil
numberLayers(coil::Coil) = trunc(Int,coil.thickness/(2*coil.wireRadius)) #Number of winding layers in the coil
totalNumberWindings(coil::Coil) = trunc(Int,numberLayers(coil)*(numberWindings(coil)-0.5))
wireLength(coil::Coil) = pi*numberLayers(coil)*(numberWindings(coil)*coil.innerRadius+sqrt(3)*coil.wireRadius*(numberWindings(coil)*(numberLayers(coil)+1)-(numberLayers(coil)+3)/2))
wireArea(coil::Coil) = coil.wireRadius^2*pi
wireVolume(coil::Coil) = wireLength(coil)*wireArea(coil.wireRadius)
wireMass(coil::Coil,proj::Projectile) = density(proj)*wireVolume(coil)
resistance(coil::Coil) = resistivityCu*wireLength(coil)/wireArea(coil)
magDomainVol(proj::Projectile) = proj.domainSize^3
coilCrossSectionalArea(coil::Coil) = (coil.innerRadius+coil.thickness) * coil.length

function magneticFieldSummation(coil::Coil, current::Current, positionFromCoil::Length)
    coilRadius = coil.innerRadius
    wireRad = coil.wireRadius
    μ0/2*current*sum((layerNumber*sqrt(3)*wireRad+coilRadius)^2/((layerNumber*sqrt(3)*wireRad+coilRadius)^2+(positionFromCoil-(2*wireRad*(rowNumber-numberWindings(coil)/2)))^2)^(3/2) for rowNumber=1:numberWindings(coil) for layerNumber=1:numberLayers(coil)) |>ustrip
end

#Reminder: The point starts in middle of the coil, then moves outward and goes through the front of the coil. When intPostion = coilLength it's at CoilFront.
function magneticFieldIntegration(coil::Coil, current::Current, integrationPositon::Length) :: BField
    coilInnerRadius = coil.innerRadius
    coilOuterRadius = coil.innerRadius + coil.thickness
    crossSectionalArea = coilCrossSectionalArea(coil)
    distToCoilBack = integrationPositon + coil.length/2
    distToCoilFront = distToCoilBack - coil.length
    constant = μ0*totalNumberWindings(coil)*current*0.25
    logarithm(coilPosition::Length,radialLength::Length) = log((coilPosition|>ustrip)^2+(radialLength|>ustrip)^2)
    arctan(coilPosition::Length,radialLength::Length) = atan((coilPosition/radialLength)|>ustrip)
    intigration(radialLength::Length) = 
        distToCoilBack  * logarithm(distToCoilBack,radialLength) - 
        distToCoilFront * logarithm(distToCoilFront,radialLength)+
        (distToCoilBack-distToCoilFront) + 2 * radialLength * 
        (arctan(distToCoilBack,radialLength) - 
        arctan(distToCoilFront,radialLength))

    return (constant * (intigration(coilOuterRadius) - intigration(coilInnerRadius)) / crossSectionalArea) |> T
end
function magneticFieldIntegration(coil::Coil, current::Current, coilPosition::Length, globalPosition::Length) :: BField
    magneticFieldIntegration(coil, current, coilPosition - globalPosition)
end


export IronProjectile, NickelProjectile, Coil, Barrel, volume, mass, density, numberWindings, numberLayers, wireLength, wireArea, wireVolume, wireMass, resistance, magDomainVol,magneticFieldSummation,magneticFieldIntegration
end #module