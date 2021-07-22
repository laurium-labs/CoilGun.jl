using CoilGun: Barrel, CoilGenerator, resistance,  Scenario, IronProjectile, ProjectilePhysical, ProjectileMagnetic, ProjectileCoilEvent, solveScenario, Projectile, InitialConditions
using CoilGun: Coil, Barrel, ProjectileCoilEvent, Voltage, ElectricalResistance, HField, Length, Velocity
using Unitful:Ω, m, cm, kg, g, A, N, Na, T, s, μ0, ϵ0, k, J, K, mol, me, q, ħ, μB, mm, inch, μm, H, V, gn
using Unitful:Length, Mass, Current, Capacitance, Charge, Force, ElectricalResistance, BField, Volume, Area, Current, HField, MagneticDipoleMoment, Density, Inductance, ustrip, Voltage, Acceleration, Time, Velocity
using ForwardDiff
using DifferentialEquations.EnsembleAnalysis: timestep_mean
#using Plots
using Test
using Plots
include("scenario_constants.jl")
include("set_up_scenario.jl")
# Curr = map(i ->CoilGun.current(coils[i], totalΩ, volts, t - eventTimes[1].entersActiveZone, magnetization, velocity, position), 1:length(coils))
# I = Curr[1]
# B = bFieldCoil(coil, I, position)
# ∇B = ∇BFieldCoil(coil, I, position)
# Magirr = Mag_irr(ip, B, Magirr, magnetization)
# dH = dHField(coils, volts, totalΩ,∇B, magnetization, position, velocity, accel, t) 
# magnetization += ∂Magnetization_∂HField(ip, B, Magirr, dH) * dH * Δt

# ∂current = ∂Current(coil, t, volts, totalΩ, position, velocity, acceleration(totalForce(ip, ∇B, velocity, magnetization), mass(ip)), magnetization)
# ∂BField_∂Current(coil, position)
# include("unitTests.jl")
# return
# println("Finished Unit tests")

tspan = (0.0,0.16)
tsteps = range(tspan[1], tspan[2], length=50)
initial = InitialConditions(magIrr, position, velocity, magnetization)
scenario = Scenario(
    ip,
    bar,
    coils,
    PCE,
    tspan,
    tsteps,
    volts,
    resistor,
    initial
)

println("Now Solving...")
sln = solveScenario(scenario)
println("Length of Velocity:\t\t",length(sln[4,:]),"\nLength of Position:\t\t", length(sln[3,:]),"\nLength of Time:\t\t\t", length(sln[2,:]),"\nLength of Magnetization:\t", length(sln[1,:]))
#figure(figsize=(8,6))
p1 = plot(sln, vars=(0,1), ylabel = "[m]"  , title = "Displacement", legend = false)
p2 = plot(sln, vars=(0,2), ylabel = "[m/s]", title = "Velocity", legend = false)
p3 = plot(sln, vars=(0,4), ylabel = "[A/m]", title = "Magnetization", legend = false)
p4 = plot(sln, vars=(0,3), ylabel = "[A/m]", title = "Irriversible Magnetization", legend = false)
display(plot(p1,p2,p3,p4, layout = (2,2)))

# dist = coilLen|>m|>ustrip
# println("Max Velocity $(sln[3,:][argmax(sln[3,:])])m/s")
# println("Point were projectile started to accerate $(dist .- sln[2,:][argmin(sln[3,:])])m.\nPoint where projectile started to decelerate $(sln[2,:][argmax(sln[3,:])] .- dist)m")
# plot(sln, layout = (2,2))