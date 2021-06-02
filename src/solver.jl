using DifferentialEquations: ODEProblem, solve

function coilTime(t::Time, eventTimes::ProjectileCoilEvent, coilNumber::Int)::Time
    return if !isnothing(eventTimes.exitsActiveZone[coilNumber]) && t - eventTimes.exitsActiveZone[coilNumber] >= 0s
        # println("exitsActiveZone works")
        t - eventTimes.exitsActiveZone[coilNumber]|>s
    elseif !isnothing(eventTimes.entersActiveZone[coilNumber]) && t - eventTimes.entersActiveZone[coilNumber] >= 0s
        # println("entersActiveZone works")
        t - eventTimes.entersActiveZone[coilNumber] |>s 
    else
        t
    end
end

struct Scenario
    proj::Projectile
    barrel::Barrel
    coils::Array{Coil,1}
    eventTimes::ProjectileCoilEvent
    endTime::Time
    voltage::Voltage
    resistor::ElectricalResistance
    initialMagIRR::HField
    initalPosition::Length
    initialVelocity::Velocity
    initialMagnetization::HField
end 

function coilProblem!(du,u,scenario,time)
    position        = view(u , 1)m
    ∂Position_∂t    = view(du, 1)
    velocity        = view(u , 2)m/s
    acceleration    = view(du, 2)
    magIrr          = view(u , 3)A/m
    ∂MagIrr_∂t      = view(du, 3)
    magnetization   = view(u , 4)A/m
    ∂Mag_∂t         = view(du, 4)
    t               = (time)s

    totalΩ = scenario.resistor + resistance(scenario.coils[1])
    foreach(eachindex(scenario.coils)) do coilIndex
        distanceFromCoil = scenario.coils[coilIndex].location - position
        reasonable = velocity/coilIndex < 10m/s
        if reasonable
            if distanceFromCoil <= scenario.coils[coilIndex].effectiveRange && isnothing(scenario.eventTimes.entersActiveZone[coilIndex])
                println("\nCoil $(coilIndex) turned on")
                scenario.eventTimes.entersActiveZone[coilIndex] = t
            end
            if distanceFromCoil <= 0.0m && isnothing(scenario.eventTimes.exitsActiveZone[coilIndex])
                println("Exited coil\t$(coilIndex)\t with velocity:\t $(velocity) at time $(t),\t acceleration: $(acceleration[1])")
                scenario.eventTimes.exitsActiveZone[coilIndex] = t
            end
        end
    end
    curr = map(i -> current(scenario.coils[i], totalΩ, scenario.voltage,coilTime(t,scenario.eventTimes,i),magnetization,velocity,position),eachindex(scenario.coils))
    H  = sum(hFieldCoil( scenario.coils[i], curr[i], position) for i in eachindex(curr))
    ∇H = sum(∇HFieldCoil(scenario.coils[i], curr[i], position) for i in eachindex(curr))
    dH = dHField(scenario.coils, scenario.voltage, totalΩ, ∇H, position, velocity, scenario.eventTimes, t)
    dHe_dH = 1 + α * ∂Magnetization_∂HField(scenario.proj, H, magnetization, magIrr, dH)

    ∂Position_∂t[:] .= velocity |> ustrip
    acceleration[:] .= totalForce(scenario.proj, ∇H, velocity, magnetization)/mass(scenario.proj) |> ustrip
    ∂MagIrr_∂t[:]   .= ∂Mag_irr_∂He(scenario.proj, H, magnetization, magIrr, dH) * dHe_dH * dH |> A/(m*s) |> ustrip
    ∂Mag_∂t[:]      .= ∂Magnetization_∂HField(scenario.proj, H, magnetization, magIrr, dH) * dH |> A/(m*s) |> ustrip
    return nothing
end

function solveScenario(scenario::Scenario) 
    u0 = ustrip.([
        scenario.initalPosition      |> m, 
        scenario.initialVelocity     |> m/s, 
        scenario.initialMagIRR       |> A/m,
        scenario.initialMagnetization|> A/m 
    ])
    tspan = (0.0, scenario.endTime |> s |> ustrip)
    problem = ODEProblem(coilProblem!, u0,tspan, scenario)
    return solve(problem, alg_hints = [:stiff])
end