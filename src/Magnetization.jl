#The paper referenced for these following equations relating to the magnetization of the projectile makes use of the Wiess mean Field theory in order to predict how the sample as a whole will react under a certain magnetic field.
δ(inc::CreatedUnits.HFieldRate) = iszero(inc) ? 1 : sign(inc)
function δM(proj::Projectile, hField::HField, magnetization::HField, mag_Irr::HField, inc::CreatedUnits.HFieldRate)::Int
    #This corrects for when the field is reversed, and the difference between the irriversible magnetization (Mag_irr) and the and the anhysteris magnetization is the reversible magnetization. This function should take the values of 1 or 0.
    Mrev = ℒ(proj, hField, magnetization) - mag_Irr
    isnan(Mrev) && (@show @__LINE__,ℒ(proj, hField, magnetization), hField, magnetization)
    # println("δM:\tMrev $(Mrev)")
    return iszero(Mrev) ? 1 : 0.5 + sign(Mrev * δ(inc))/2
end

function ℒ(proj::Projectile, hField::HField, magnetization::HField)::HField
    #langevin funciton that represents the anhystesis bulk magnetization for a given material. It can be imagined as a sigmoid shape on a M-H graph.
    constant = k*roomTemp/magMomentPerDomain |>A/m  |> ustrip
    effectiveHField = hField + proj.magnetic.interdomainCoupling * magnetization |> A/m |> ustrip
    taylorApproxℒ(x) = x/(3*constant) - x^3/(45*constant^3)
    # effectiveHField == 0.0 && (@show taylorApproxℒ(effectiveHField))
    # isnan((coth(effectiveHField/constant) - constant/effectiveHField) * proj.magnetic.saturationMagnetization) && (@show hField, magnetization)
    return iszero(effectiveHField) ? 0.0A/m : (coth(effectiveHField/constant) - constant/effectiveHField) * proj.magnetic.saturationMagnetization
end

function ∂ℒ(proj::Projectile, hField::HField, magnetization::HField)::Float64
    #The first order derivative (with respect to the BField) of the ℒ function
    constant = k*roomTemp/magMomentPerDomain |>A/m|>ustrip             #Constant
    effectiveHField = hField + proj.magnetic.interdomainCoupling * magnetization |>A/m|>ustrip  #Variable
    # ∂taylorApproxℒ(x) = 1/(3*constant) - x^2/(15*constant^3)
    langevin(x) = coth(x/constant) - constant/x
    return iszero(effectiveHField) ? 0.0 : ForwardDiff.derivative(langevin,effectiveHField)*m/A * proj.magnetic.saturationMagnetization
end

function mag_Irr(proj::Projectile, hField::HField, magnetization::HField)::HField
    #This calculates the bulk irriversible magnetization inside the projectile.
    return (magnetization - proj.magnetic.reversibility * ℒ(proj,hField,magnetization))/(1-proj.magnetic.reversibility)
end

function ∂Mag_irr_∂He(proj::Projectile, hField::HField, magnetization::HField, mag_Irr::HField, dH::CreatedUnits.HFieldRate)::Float64
    # println("∂Mag_irr_∂He:\tℒ $(ℒ(proj,hField,mag_Irr)),\tmag_Irr $(mag_Irr)")
    return (ℒ(proj,hField,magnetization) - mag_Irr)/(domainPinningFactor * δ(dH))
end

function ∂Magnetization_∂HField(proj::Projectile, hField::HField, magnetization::HField, mag_Irr::HField, dH::CreatedUnits.HFieldRate)::Float64
    #Change in the objects magnetization due to an external B-Field.
    M_rev = δM(proj,hField,magnetization, mag_Irr,dH) * (ℒ(proj,hField,magnetization) - magnetization)
    numerator = M_rev + proj.magnetic.reversibility * ∂ℒ(proj, hField, magnetization) * (domainPinningFactor * δ(dH))
    # println("∂Magnetization_∂HField:\tδ: $(δ(∂H)),\t∂ℒ:$(∂ℒ(proj, hField, Mag_irr)),\tΔM_irr:$(ΔM_irr),\tnumerator: $(numerator)")
    denominator = (domainPinningFactor * δ(dH)) - proj.magnetic.interdomainCoupling * numerator
    return numerator/denominator
end