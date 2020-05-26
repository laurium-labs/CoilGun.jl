using DifferentialEquations

struct Scenario
    proj::Projectile
    barrel::Barrel
    coils::Array{Coil,1}
    endTime::Time
    voltage::Voltage
    resistor::ElectricalResistance
    initialMagIRR::HField
    initalPosition::Length
    initialVelocity::Velocity
    initialMagnetization::HField
end 

function coilProblem!(du,u,scenario,time )
    magnetization = (u[1])A/m
    ∂Mag_∂t = view(du, 1)
    position = (u[2])m
    ∂Position_∂t = view(du, 2)
    velocity = (u[3])m/s
    acceleration = view(du, 3)
    magIrr = (u[4])A/m
    ∂MagIrr_∂t = view(du,4 )
    t = (time)s

    totalΩ = scenario.resistor + resistance(scenario.coils[1])
    I = map(coil -> current(coil, totalΩ, scenario.voltage, t, magnetization, velocity, position), scenario.coils)

    B = sum(map(i -> bFieldCoil(scenario.coils[i], I[i], position), 1:length(I)))
    ∇B = sum(map(i -> ∇BFieldCoil(scenario.coils[i], I[i], position), 1:length(I)))
    force = totalForce(scenario.proj, ∇B, velocity[1], magnetization[1])
    accel = (force/mass(scenario.proj)) |> m/(s^2)
    acceleration[1] = (accel) |> ustrip
    ∂Position_∂t[1] = velocity |> m/s |> ustrip
    dH = ∂HField(scenario.coils, I, scenario.voltage, totalΩ,∇B, magnetization, position, velocity, accel, t) 
    ∂Mag_∂t[1] = ∂Magnetization_∂HField(scenario.proj, B, magIrr, dH) * dH |> A/(m*s) |> ustrip
    ∂MagIrr_∂t[1] = ∂Mag_irr_∂He(scenario.proj,δ(dH), δM(scenario.proj, B, magIrr, dH), ℒ(scenario.proj, B, magIrr), magIrr) * dH |> A/(m*s) |> ustrip
    nothing
end

function solveScenario(scenario::Scenario) 
    u0 = [scenario.initialMagnetization |> A/m |> ustrip, scenario.initalPosition |> m |> ustrip, scenario.initialVelocity |> m/s |> ustrip, scenario.initialMagIRR |> A/m |> ustrip]
    tspan = (0.0, scenario.endTime |> s |> ustrip)
    problem = ODEProblem(coilProblem!, u0,tspan, scenario)
    return solve(problem)
end