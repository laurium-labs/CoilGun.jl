using CoilGun
using Unitful:Ω, m, cm, kg, g, A, N, Na,T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm
using Unitful:ustrip, Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density
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
const roomTemp = 293K                                       #Standard room Tempearture
const domainPinningFactor = 150A/m                          #This is the domain pinning factor for Iron (transformer)
const domainMagnetization = 0.2 * numberAtomsperDomainFe*bohrMagnetonPerAtomFe |> A/m #Magnetization of the domain
const magMomentPerDomain = domainMagnetization*domainSizeFe^3    #This dipole magnetic moment doesn't take hysteresis/pinning into effect
const saturationMagnetizationPerKgFe = 217.6A/(m*kg)             #Saturation magnetizaiton of pure Iron per unit mass.


#Projectile Specifications
projrad = 3.5mm |> m
projlength = 1inch |> m
magdomiansize = 26.5μm |> m   #From literature, it is actually 26.5 nm, this value is not being used because I would have no memory left.
magstrngth = 1T
saturationMagnetization = saturationMagnetizationPerKgFe * magdomiansize^3 * densityFe

#Barrel Specifications
bthickness = 1mm |> m #Barrel thickness
blength = 0.5m |> m #Length of barrel


#Coil Specifications
innerRadius = projrad+bthickness      #The Inner diameter of the Coil needs to be the Outer diameter of the barrel
cthickness = 1inch |> m      #The difference in the Inner diameter and Outer diameter of the Coil
coilLen = projlength    #The length of the Coil should be the exact length of the projectile
coilHght = 2.3e-2m |> m      #Distance from inner to outer diameter of the Coil
wirerad = 1.6mm |> m         #The radius of 14-guage wire including insulation
position = projlength

I = 1A #Current flowing through the wire
stepSize = 1_000

phys    = ProjectilePhysical(projrad,projlength,densityFe)
mag     = ProjectileMagnetic(domainSizeFe,magstrngth,α,magMomentPerDomain,domainMagnetization,saturationMagnetization,generateMagneticDomians(phys,domainSizeFe,μ0 * saturationMagnetization),0T)
ip      = IronProjectile(phys,mag,position)
bar     = Barrel(ip.physical.radius,bthickness,blength)
coil    = Coil(bar.innerRadius+bar.thickness,cthickness,ip.physical.length,wirerad)

"""
When calculating the magnetic field given off from the coil, remember to have the step size similar to the domain size in the projectile.
"""
radii = map(step -> step * (coil.length*3/2) / stepSize , 1:stepSize)
# bFieldSummation = map(radii) do positionFromCoil
#     return magneticFieldSummation(coil, I, positionFromCoil) |>ustrip
# end

bFieldIntegration = map(radii) do positionFromCoil
    return magneticFieldIntegration(coil, I, positionFromCoil)
end

magmin = 0A/m
magmax = 1A/m
bf = 0.3T
pbf = 0.2T
mem = [0T, 1T, 0.2T]
mem = Array{BField,1}(mem)
closingFunction(magmin,magmax,ip,bf,pbf)
effectiveMagnetization(ip, bf, mem)

bfg1 = 2T/m
testbfg = hcat([0T/m : 0.001T/m : bfg1]...)
bfg = BFieldGradient(testbfg)
dipoleCoilForce([12,14], ip, coil, bfg1)
projectileCoilTotalForce(coil, ip, bfg)