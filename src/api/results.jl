
# function solved_default(default_scenario::Scenario)
#     return solveScenario(default_scenario)
# end

function extract_results(sln)
    # println("Length of Velocity:\t\t",length(sln[4,:]),"\nLength of Position:\t\t", length(sln[3,:]),"\nLength of Time:\t\t\t", length(sln[2,:]),"\nLength of Magnetization:\t", length(sln[1,:]))
    # sln=solved_default(default)
    time=sln[0,:]
    velocity = sln[3,:]
    magentization= sln[1,:]
    irreversibleMagentization=sln[4,:]
    displacement = sln[2,0]

    #body = Dict("time"=>stripped_time_core, "temperature" => stripped_temperature_core)
    #result_dict = Dict("time" => time, "displacement" => displacement, "velocity" => velocity, "magentization" => magentization, "irreversibleMagentization"=>irreversibleMagentization)
    return string(time, velocity, magentization, irreversibleMagentization, displacement)
end