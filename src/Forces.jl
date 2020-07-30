function dipoleCoilForce(proj::Projectile, ∇HField::CreatedUnits.HFieldGrad, magnetization::HField)::Force #Update: Integrate mag over volume
    #This force function assumes a number of things:
    #1) The magnetization of the projectile is constant throughout the material. This means that all of the magnetic domains are consistent throughout the material.
    #2) The magnetic field experienced at the center of the projectile is the average magnetic field experienced by the projectile.
    #3) The magnetization of the material can be approximated as a magnetic dipole (loop would be more accurate, but this is easier).
    magneticDipoleMoment = magnetization * volume(proj)
    return magneticDipoleMoment * μ0 * ∇HField|>N
end
function frictionForce(proj::Projectile, velocity::Velocity, dipoleCoilForce::Force)::Force #Update (check inside)
    #This describes the resistive force due to friction (both static and kinetic). I know very little about the relationship between kinetic friction and velocity, but I highly doubt it's a constant relationship. More research will have to be done.
    frictionCoefficient = abs(velocity) > 1e-5m/s ? kineticFrictionCoefficientFe : staticFrictionCoefficientFe
    normalForce = mass(proj) * gn
    signCorrection = velocity == 0 ? 1 : -sign(velocity)
    if abs(normalForce) >= abs(dipoleCoilForce) && abs(velocity) < 1e-5m/s
        return -dipoleCoilForce|>N
    else
        return frictionCoefficient * normalForce * signCorrection|>N
    end
end
function airResistance(proj::Projectile, velocity::Velocity)::Force #Update to model rod through tube
    #This funciton is used to calculate the air resistance on the projectile. This current function is overly simplified and will need to be changed later for a more accurate function.
    return 6 * pi * dynamicViscosityAir * proj.physical.radius * velocity|>N
end
function totalForce(proj::Projectile, ∇HField::CreatedUnits.HFieldGrad, velocity::Velocity, magnetization::HField)::Force
    dipoleForce = dipoleCoilForce(proj,∇HField, magnetization)
    frictionalForces = frictionForce(proj, velocity, dipoleForce) + airResistance(proj, velocity)
    return dipoleForce + frictionalForces |> N
end