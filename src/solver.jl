using DifferentialEquations

function coilTime(t::Time, eventTimes::ProjectileCoilEvent, ind::Int)::Time
    if !isnothing(eventTimes.exitsActiveZone[ind]) && t - eventTimes.exitsActiveZone[ind] >= 0s
        # println("exitsActiveZone works")
        return t - eventTimes.exitsActiveZone[ind]|>s
    elseif !isnothing(eventTimes.entersActiveZone[ind]) && t - eventTimes.entersActiveZone[ind] >=0s
        # println("entersActiveZone works")
        return t - eventTimes.entersActiveZone[ind] |>s 
    else
        return t
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

function coilProblem!(du,u,scenario,time )
    magnetization   = (u[1])A/m
    ∂Mag_∂t         = view(du, 1)
    position        = (u[2])m
    ∂Position_∂t    = view(du, 2)
    velocity        = (u[3])m/s
    acceleration    = view(du, 3)
    magIrr          = (u[4])A/m
    ∂MagIrr_∂t      = view(du,4 )
    t               = (time)s

    totalΩ = scenario.resistor + resistance(scenario.coils[1])
    foreach(1:length(scenario.coils)) do coilInd
        distanceFromCoil = scenario.coils[coilInd].location - position
        if distanceFromCoil <= scenario.coils[coilInd].coilOnRange && isnothing(scenario.eventTimes.entersActiveZone[coilInd])
            println("\nEntered coil\t$(coilInd)")
            scenario.eventTimes.entersActiveZone[coilInd] = t
        end
        if distanceFromCoil <= 0.0m && isnothing(scenario.eventTimes.exitsActiveZone[coilInd])
            println("Exited coil\t$(coilInd)\t with velocity:\t $(velocity),\t and coil efficiency of: $(velocity/coilInd)")
            scenario.eventTimes.exitsActiveZone[coilInd] = t
        end
    end
    curr = map(i -> current(scenario.coils[i], totalΩ, scenario.voltage,coilTime(t,scenario.eventTimes,i),magnetization,velocity,position),1:length(scenario.coils))
    H  = sum(map(i ->  hFieldCoil(scenario.coils[i], curr[i], position), 1:length(curr)))
    ∇H = sum(map(i -> ∇HFieldCoil(scenario.coils[i], curr[i], position), 1:length(curr)))
    # println("∇H: $(H), v: $(velocity)")
    force = totalForce(scenario.proj, ∇H, velocity[1], magnetization[1])
    accel = (force/mass(scenario.proj)) |> m/(s^2)
    acceleration[1] = (accel) |> ustrip
    ∂Position_∂t[1] = velocity |> m/s |> ustrip
    dH = dHField(scenario.coils, scenario.voltage, totalΩ, ∇H, position, velocity, scenario.eventTimes, t)
    # println("force $(force),\tmagnetization $(magnetization[1]),\tdH $(dH)")
    ∂MagIrr_∂t[1] = ∂Mag_irr_∂He(scenario.proj,H, magIrr, dH) * dH |> A/(m*s) |> ustrip
    ∂Mag_∂t[1] = ∂Magnetization_∂HField(scenario.proj, H, magIrr, dH) * dH |> A/(m*s) |> ustrip
    # println("∂Mag_∂t:$(∂Mag_∂t[1]),\t∂Magnetization_∂HField:$(∂Magnetization_∂HField(scenario.proj, H, magIrr, dH)),\tdH:$(dH)")
    nothing
end

function solveScenario(scenario::Scenario) 
    u0 = [scenario.initialMagnetization |> A/m |> ustrip, scenario.initalPosition |> m |> ustrip, scenario.initialVelocity |> m/s |> ustrip, scenario.initialMagIRR |> A/m |> ustrip]
    tspan = (0.0, scenario.endTime |> s |> ustrip)
    problem = ODEProblem(coilProblem!, u0,tspan, scenario)
    return solve(problem, alg_hints = [:stiff])
end