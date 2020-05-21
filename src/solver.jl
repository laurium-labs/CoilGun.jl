using DifferentialEquations

struct Scenario
    proj::Projectile
    barrel::Barrel
    coil::Coil
    endTime::Time
    voltage::Voltage
    resistor::ElectricalResistance
    initialCurrent::Current
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

    totalΩ = scenario.resistor + resistance(scenario.coil)
    current = current(scenario.proj, scenario.coil, totalΩ, scenario.voltage, t, magnetization, velocity, position)
    B = simpleBField(scenario.coil, current, position)
    ∇B = bFieldGradient(scenario.coil, I, position)
    Magirr = Mag_irr(scenario.proj, B, Magirr, magnetization)
    dH = ∂HField(scenario.coil, current, volts, totalΩ,∇B, magnetization, position, velocity, accel, t) 
    # (∇B * ip.velocity + simpleBField(coil, I - prevI, ip.position)/t) * t / μ0
    magnetization += ∂Magnetization_∂HField(ip, B, Magirr, dH) * dH * Δt
    
end

function solveScenario(scenario::Scenario) 
    u0 = [ scenario.initialMagnetization, scenario.initalPosition, scenario.initialVelocity]
end