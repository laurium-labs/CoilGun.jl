
#=
    This file is strictly for testing the api folder with CoilGun
=#



using CoilGun: dictionary_api, get_default_scenario_json, Scenario
const barrel_Length = "CCoilGun.Scenario(CoilGun.IronProjectile(CoilGun.ProjectilePhysical(0.0035 m, 56000//1000 m, 7750 kg m⁻³), CoilGun.ProjectileMagnetic(2.65e-6 m, 0.00134, 1.61e6 A m⁻¹, 0.373)), CoilGun.Barrel(1//1000 m, 1//100 m, 0.5 m), CoilGun.Coil[CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.01115 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.0223 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.03345 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.0446 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.05575 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.0669 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.07805000000000001 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.0892 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.10035 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.1115 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.12265000000000001 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.1338 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.14495 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.15610000000000002 m, 0.0223 m), CoilGun.Coil(0.0035 m, 0.0135 m, 0.01115 m, 0.0016 m, 0.16725 m, 0.0223 m)], CoilGun.ProjectileCoilEvent(Union{Nothing, Union{Unitful.Quantity{T,𝐓,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐓,U}} where S where L} where U where T}[nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing], Union{Nothing, Union{Unitful.Quantity{T,𝐓,U}, Unitful.Level{L,S,Unitful.Quantity{T,𝐓,U}} where S where L} where U where T}[nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing]), 1 s, 15 V, 10 Ω, 1 A m⁻¹, 0 m, 0.1 m s⁻¹, 0 A m⁻¹)"
const coil_Length = "coilLength" => 100
const barrel_thickness = "barrelThickness" => 100
const coil_wire_radius = "coilWireRadius" => 100
# dictionary_api(Dict("initialMagIRR"=>10, 
# "proj.magnetic.saturationMagnetization"=>1610000,
# "proj.physical.radius"=>0.004,
# "endTime"=>4,
# "initialMagnetization"=>0,
# "barrel.thickness"=>0.399,
# "barrel.innerRadius"=>0.001,
# "resistor"=>10,
# "proj.physical.density"=>7751,
# "voltage"=>15,
# "initalPosition"=>0,
# "proj.physical.length"=>0.054,
# "barrel.length"=>45,
# "proj.magnetic.domainSize"=>0.00000265,
# "initialVelocity"=>10))
 println(get_default_scenario_json())