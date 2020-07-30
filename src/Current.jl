#Equations relating to the calculation of current
function selfInductance(coil::Coil)::Inductance #Fix for new coil Equation
    #This is a simplified version for the self inductance of the coil. It is not taking into consideration the thickness (different layers) of the coil, or the varying magnetic fields that pass through each loop. This is just for approximation only.
    iR = coil.innerRadius
    oR = coil.outerRadius
    length = coil.length
    B_I = pi*μ0*totalNumberWindings(coil)/(3*coilCrossSectionalArea(coil)^2)
    var = ((oR^2+length^2)^(3/2)+(iR^2+length^2)^(3/2)-(oR^3-iR^3))*(oR^2-iR^2)|> m^5
    return B_I*var
end
function projectileInducedVoltage(coil::Coil, magnetization::HField, velocity::Velocity, position::Length)::Voltage #Update:Includes face of coil, not depth
    position = position - coil.location
    area = pi * (coil.outerRadius^2 - coil.innerRadius^2)
    effectiveRadius = area/(coil.outerRadius - coil.innerRadius)
    ∂AreaRatio_∂t = effectiveRadius * velocity * position/(position^2 + effectiveRadius^2)^(3/2)
    constant = μ0 * magnetization * totalNumberWindings(coil) * area #includes #of windings because the BField is traveling through the windings.
    return constant * ∂AreaRatio_∂t
end
function ∂projectileInducedVoltage(coil::Coil, position::Length, velocity::Velocity, acceleration::Acceleration, magnetization::HField)#Update:Includes face of coil, not depth
    #This funciton describes the change in the induced voltage per change in time
    position = position - coil.location
    area = pi * (coil.outerRadius^2 - coil.innerRadius^2)
    effectiveRadius = area/(coil.outerRadius - coil.innerRadius)
    ∂∂AreaRatio_∂tt = effectiveRadius * (position*acceleration/(position^2+effectiveRadius^2)^(3/2) + effectiveRadius^2*velocity^2/(position^2+effectiveRadius^2)^(5/2)) |> s^-2
    constant = μ0 * magnetization * totalNumberWindings(coil) * area
    return constant * ∂∂AreaRatio_∂tt |> V/area
end
function couplingFactor(coil::Coil)::Float64 #May need to be changed. There's another equation with number of turns
    #This function calculates the ratio of magnetic field lines passing through the generating coil, and an adjacent coil.
    a = coil.length/2
    α1 = 4 * a
    β1 = 2 * a
    β2 = 0m
    block(var::Length) = (coil.outerRadius^2 + var^2)^(3/2) - (coil.innerRadius^2 + var^2)^(3/2)
    numerator = block(α1) - 2*block(β1) + block(β2)
    denominator = 2 * (block(β1)-block(β2))
    return numerator/denominator
end

#TODO: Create a function that describes how a capacitor will supply voltage to a coil
#Functions for current
function coilCurrent(coil::Coil, position::Length, time::Time, maxVoltage::Voltage, characteristicTime::Time, resistance::ElectricalResistance)::Current #Fix: Turn coil off when proj. is > 2 coilOnRange away from coil, negative time after coil reverses current
    #The entry fields are intentionally left without specification due to its derivative being taken.
    distFromCoil = coil.location - position
    if (0m <= distFromCoil) && (distFromCoil <= coil.coilOnRange)
        return maxVoltage * (1-exp(-time / characteristicTime)) / resistance
    elseif distFromCoil < 0m
        return maxVoltage * (2*exp(-time / characteristicTime)-1) / resistance
    else
        return 0A
    end
end
function current(coil::Coil, totalΩ::ElectricalResistance, initialVoltage::Voltage, time::Time, magnetization::HField, velocity::Velocity, position::Length)::Current
    #This function calculates the current that is traveling through a coil. This is not taking operational amplifiers into consideration.
    𝓀 = couplingFactor(coil)
    τ = selfInductance(coil)*𝓀^2/totalΩ #Characteristic Time (current reaches steady state value at 5*τ)
    projectileInducedCurrent = projectileInducedVoltage(coil, magnetization, velocity, position)/totalΩ
    return coilCurrent(coil, position|>m, time,initialVoltage,τ,totalΩ)|>A #+ projectileInducedCurrent |> A
end
function ∂Current(coil::Coil, time::Time, initialVoltage::Voltage, totalΩ::ElectricalResistance, position::Length)
    #This function describes how the current through the coil changes with the change in time.
    τ = selfInductance(coil)*couplingFactor(coil)^2/totalΩ #Characteristic Time (Current reaches it's steady state value at 5*τ)
    distFromCoil = coil.location - position
    rateOfChange = if (0m <= distFromCoil) && (distFromCoil <= coil.coilOnRange)
        exp(-time / τ)/τ
    elseif distFromCoil < 0m
        -2 * exp(-time / τ)/τ
    else
        0/s
    end
    # println("∂Current: ∂CoilVoltageRate $(∂CoilVoltageRate)")
    return initialVoltage * rateOfChange/totalΩ|>A/s# + ∂projectileInducedVoltage(coil,position,velocity,acceleration,magnetization))/totalΩ |> A/s
end