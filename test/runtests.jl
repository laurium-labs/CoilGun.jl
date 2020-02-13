using CoilGun
using Unitful:Ω,m,mm,cm,μm,inch,kg,A,N,ustrip,T,s,μ0,ϵ0,Length,Mass,Current,Capacitance,Charge,Force,ElectricalResistance,BField
using Test

#Projectile Specifications
projrad = 3.5mm
projlength = 1inch
magdomiansize = 26.5μm

#Barrel Specifications
bthickness = 1mm #Barrel thickness
blength = 0.5m #Length of barrel

#Coil Specifications
innerRadius = projrad+bthickness      #The ID of the Coil needs to be the OD of the barrel
cthickness = 1inch      #The difference in the ID and OD of the Coil
coilLen = projlength    #The length of the Coil should be the exact length of the projectile
coilHght = 2.3e-2m      #Distance from inner to outer diameter of the Coil
wirerad = 1.6mm         #The radius of 14guage wire including insulation

I = 1A #Current flowing through the wire
stepSize = 1_000

ip = IronProjectile(projrad,projlength,magdomiansize)
bar = Barrel(ip.radius,bthickness,blength)
coil = Coil(bar.innerRadius+bar.thickness,cthickness,ip.length,wirerad)

radii = map(step -> step * (coil.length*3/2) / stepSize , 1:stepSize)
bFieldSummation = map(radii) do positionFromCoil
    return magneticFieldSummation(coil, I, positionFromCoil)
end

bFieldIntegration = map(radii) do positionFromCoil
    return magneticFieldIntegration(coil, I, positionFromCoil)
end
    

println(bFieldIntegration[1],"\n",bFieldIntegration[end])