using CoilGun: dictionary_api
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn
using Unitful:Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Acceleration, Time, Velocity
barrelLength = "barrelLength"
barrelThickness = "barellThickness"
wireRadius = "wireRadius"
coilLength = "coilLength"

treeDict = Dict(barrelLength => "4", wireRadius => "4", coilLength => "4", barrelThickness => "4" )
println(treeDict)
dictionary_api(treeDict)