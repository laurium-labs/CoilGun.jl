using DifferentialEquations

struct Scenario
    proj::Projectile
    barrel::Barrel
    coil::Coil
    endTime::Time
    voltage::Voltage
    resistor::ElectricalResistance
    initialMagIRR::HField
    initalPosition::Length
    initialVelocity::Velocity
    initialMagnetization::HField
end 

function coilProblem!(du,u,scenario,t )
    magnetization = view(u, 1)
    ∂Mag_∂t = view(du, 1)
    position = view(u, 2)
    ∂Position_∂t = view(du, 2)
    velocity = view(u, 3)
    acceleration = view(du, 3)
    magIrr = view(u, 4)
    ∂MagIrr_∂t = view(du,4 )
    B = simpleBField(scenario.coil, current, position)
    ∇B = bFieldGradient(scenario.coil, current, position)

    acceleration[1] = totalForce(scenario.proj, ∇B, velocity[1], magnetization[1])
    ∂Position_∂t[1] = velocity[1]
    totalΩ = scenario.resistor + resistance(scenario.coil)
    current = current(scenario.proj, scenario.coil, totalΩ, scenario.voltage, t, magnetization, velocity, position)
    dH = ∂HField(scenario.coil, current, scenario.voltage, totalΩ,∇B, magnetization, position, velocity, acceleration[1], t) 
    ∂Mag_∂t[1] = ∂Magnetization_∂HField(scenario.proj, B, magIrr[1], dH) * dH 
    ∂MagIrr_∂t[1] = ∂Mag_irr_∂H(scenario.proj, δM(scenario.proj, B, magIrr[1], dH), ℒ(scenario.proj, B, magIrr[1]), magIrr[1]) * dH
    nothing
end

function solveScenario(scenario::Scenario) 
    u0 = [ scenario.initialMagnetization, scenario.initalPosition, scenario.initialVelocity, scenario.initialMagIRR]
    tspan = (0.0s, scenario.endTime)
    problem = ODEProblem(coilProblem!, u0,tspan)
    return solve(problem)
end