using DiffEqFlux, DifferentialEquations, Plots, GalacticOptim

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
u0 = ustrip.([
        scenario.initial.initalPosition      |> m, 
        scenario.initial.initialVelocity     |> m/s, 
        scenario.initial.initialMagIRR       |> A/m,
        scenario.initial.initialMagnetization|> A/m 
    ])

println("Now Solving initial...")
# ode_data = Array(solveScenario(scenario))
println(typeof(ode_data))

# dudt2 = FastChain((x, p) -> x.^3,
#                   FastDense(4, 50, tanh),
#                   FastDense(50, 4))

# prob_neuralode = NeuralODE(dudt2, tspan, Tsit5(), saveat = tsteps)

# dudt_no_position = FastChain((x, p) -> x.^3,
#                       FastDense(5, 50, tanh),
#                       FastDense(50, 3))
position_velocity= onehot(1, [0,1,0,0])
dudt_flux = Parallel(vcat,Dense([0 1 0 0], false, identity), Chain(Dense(4, 30, tanh), Dense(30,3)))
dudt_final3 = Chain(Dense(4, 30, tanh), Dense(30,3))
# ps = Flux.params(dudt_flux[2:2])

function neural_coil!(du, u, p, t)
    du[1] = u[2]
    du[2:4] .= dudt_final3(u, p)
end
p = randn(4*30+30+30*3+3)

prob_coilode = ODEProblem(neural_coil!, u0, tspan,p)

prob_neuralode_flux = NeuralODE(dudt_flux, tspan, Tsit5(), saveat = tsteps)
function predict_n_ode()
        _prob = remake(prob,p=p)

    Array(solve(_prob,Tsit5(),u0=z0,p=p,callback=cb,saveat=t,sensealg=ReverseDiffAdjoint()))[1:2,:]
    #Array(solve(prob,Tsit5(),u0=z0,p=p,saveat=t))[1:2,:]
end

function predict_neuralode()
    # println(length(p))
    _prob = remake(prob_coilode,p=p)

    Array(solve(_prob, Tsit5(), saveat = tsteps))
    # _prob = remake(prob_neuralode_flux,p=p)

    # _prob(u0, p)
#   prob_coilode(u0, p))
end

function loss_neuralode()
    pred = predict_neuralode()
    loss = sum(abs2, ode_data .- pred)
    return loss, pred
end

callback = function (p, l, pred; doplot = true)
  display(l)
  # plot current prediction against data
  plt = scatter(tsteps, ode_data[1,:], label = "data")
  scatter!(plt, tsteps, pred[1,:], label = "prediction")
  if doplot
    display(plot(plt))
  end
  return false
end
data = Iterators.repeated((), 200)

Flux.train!(loss_neuralode, p,data, ADAM(0.05), 
                                          cb = callback)
