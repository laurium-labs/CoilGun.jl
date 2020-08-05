#Reminder: The point starts in middle of the coil, then moves outward and goes through the front of the coil. When intPostion = coilLength it's at CoilFront.
function ∂HField_∂Current(coil::Coil, position::Length)
    #This function describes how the BField changes with respect to the change in current
    return hFieldCoil(coil, 1A, position)/1A
end
function hFieldCoil(coil::Coil, current::Current, position::Length)::HField
    constant = totalNumberWindings(coil)/coilCrossSectionalArea(coil)
    logarithm(pos::Length)::Length = pos * log((sqrt(pos^2+coil.outerRadius^2)+coil.outerRadius)/(sqrt(pos^2+coil.innerRadius^2)+coil.innerRadius))
    farEdgeofCoil = (coil.location + coil.length/2) - position
    closeEdgeofCoil = farEdgeofCoil - coil.length
    return constant * current * (logarithm(farEdgeofCoil) - logarithm(closeEdgeofCoil)) |> A/m
end

function effectiveCoilRange(coil::Coil, current::Current, position::Length)
    
end

function ∇HFieldCoil(coil::Coil, current::Current, position::Length)::CreatedUnits.HFieldGrad
    constant = current*totalNumberWindings(coil)/coilCrossSectionalArea(coil)
    outerRadius = coil.outerRadius |> m |>ustrip
    innerRadius = coil.innerRadius |>m |> ustrip
    logarithm(pos) = pos * log((sqrt(pos^2+outerRadius^2)+outerRadius)/(sqrt(pos^2+innerRadius^2)+innerRadius))
    ∇logarithm(pos::Length)::Float64 = -ForwardDiff.derivative(logarithm,pos|>m|>ustrip)
    farEdgeofCoil = coil.location - position + coil.length/2
    closeEdgeofCoil = farEdgeofCoil - coil.length
    # println("∇HFieldCoil: Current $(current),\tposition $(position),\t∇faredge $(∇logarithm(farEdgeofCoil))")
    return constant * (∇logarithm(farEdgeofCoil) - ∇logarithm(closeEdgeofCoil)) |> A/m^2
end

function dHField(coils::Array{Coil,1}, voltage::Voltage, totalΩ::ElectricalResistance, ∇H::CreatedUnits.HFieldGrad, position::Length, velocity::Velocity, eventTimes::ProjectileCoilEvent, time::Time)::CreatedUnits.HFieldRate
    #This function calculates the change in the HField due to the change in position and the change in current
    # println("dHField:\t∇H*v:$(∇H*velocity),\t∂H_∂C:$(sum(map(coil -> ∂HField_∂Current(coil,position)*∂Current(coil,time,voltage,totalΩ,position), coils)))")
    return ∇H*velocity+sum(map(i -> ∂HField_∂Current(coils[i],position)*∂Current(coils[i],coilTime(time,eventTimes,i),voltage,totalΩ,position), 1:length(coils)))|>A/m/s
end