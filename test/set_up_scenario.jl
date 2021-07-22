

#Projectile Specifications
#Physical
projrad = 3.5mm |> m
projlength = 1.0inch |> m
position = 0m   
velocity = 0.1m/s
accel = 0m/s^2
#Magnetic
saturationMagnetization = 1.61e6A/m
reversibility = 0.373
magnetization = 1e-1A/m
magIrr = 1.0A/m
println("Initial Parameters:\n\tposition:\t\t",position,
    "\n\tvelocity:\t\t", velocity,
    "\n\tacceleration:\t\t", accel, 
    "\n\tmagnetization:\t\t", magnetization)

#Coil Specifications
bthickness = 1mm |> m #Barrel thickness

innerRadius = projrad+bthickness        #The Inner diameter of the Coil needs to be the Outer diameter of the barrel
cthickness = 1inch |> m                 #The difference in the Inner diameter and Outer diameter of the Coil
coilLen = 1.0inch |> m                  #The length of the Coil should be the exact length of the projectile
coilHght = 2.3e-2m |> m                 #Distance from inner to outer diameter of the Coil
wirerad = 1.6mm |> m                    #The radius of 14-guage wire including insulation
resistor = 10Ω
volts = 15V
numberOfCoils = 10

#Barrel Specifications
blength = numberOfCoils * coilLen |> m #Length of barrel
bRadius = projrad

phys    = ProjectilePhysical(projrad,
                            projlength,
                            densityFe)
mag     = ProjectileMagnetic(domainSizeFe,
                            α,
                            saturationMagnetization,
                            reversibility)
ip      = IronProjectile(phys,mag)
bar     = Barrel(bRadius,bthickness,blength)
coils = CoilGenerator(numberOfCoils, innerRadius, innerRadius+cthickness, coilLen, wirerad)
PCE = ProjectileCoilEvent()
PCE.entersActiveZone = [nothing for _ in coils]
PCE.exitsActiveZone = [nothing for _ in coils]


# Δt=0.1s
# t = 0s
# Magirr = 0A/m
# coil = coils[1]
totalΩ = resistor + resistance(coils[1])