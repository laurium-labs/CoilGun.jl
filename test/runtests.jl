using CoilGun
using Unitful:Ω, m, cm, mm, kg, g, A, N, ustrip, T, s, μ0, ϵ0, Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, k, J, K, inch, μm
using Test
using Plots


#Projectile Specifications
projrad = 3.5mm |> m
projlength = 1inch |> m
magdomiansize = 26.5μm |> m   #From literature, it is actually 26.5 nm, this value is not being used because I would have no memory left.

#Barrel Specifications
bthickness = 1mm |> m #Barrel thickness^
blength = 0.5m |> m #Length of barrel

#Coil Specifications
innerRadius = projrad+bthickness      #The ID of the Coil needs to be the OD of the barrel
cthickness = 1inch |> m      #The difference in the ID and OD of the Coil
coilLen = projlength    #The length of the Coil should be the exact length of the projectile
coilHght = 2.3e-2m |> m      #Distance from inner to outer diameter of the Coil
wirerad = 1.6mm |> m         #The radius of 14guage wire including insulation

I = 1A #Current flowing through the wire
stepSize = 1_000

ip = IronProjectile(projrad,projlength,magdomiansize)
bar = Barrel(ip.radius,bthickness,blength)
coil = Coil(bar.innerRadius+bar.thickness,cthickness,ip.length,wirerad)


"""
When calculating the magnetic field given off from the coil, remember to have the step size similar to the domain size in the projectile.
"""
radii = map(step -> step * (coil.length*3/2) / stepSize , 1:stepSize)
bFieldSummation = map(radii) do positionFromCoil
    return magneticFieldSummation(coil, I, positionFromCoil) |>ustrip
end

bFieldIntegration = map(radii) do positionFromCoil
    return magneticFieldIntegration(coil, I, positionFromCoil)
end
