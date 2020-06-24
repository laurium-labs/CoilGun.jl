
# function solved_default(default_scenario::Scenario)
#     return solveScenario(default_scenario)
# end
function extract_results(sln)
    # println("Length of Velocity:\t\t",length(sln[4,:]),"\nLength of Position:\t\t", length(sln[3,:]),"\nLength of Time:\t\t\t", length(sln[2,:]),"\nLength of Magnetization:\t", length(sln[1,:]))
    # sln=solved_default(default)
    velocity = sln[3,:]
    magnetization= sln[1,:]
    irreversibleMagentization=sln[4,:]
    # time=sln[0,:]
     displacement = sln[2,:]

    #body = Dict("time"=>stripped_time_core, "temperature" => stripped_temperature_core)
    result_dict = Dict("velocity" => velocity, "magnetization" => magnetization, "irreversibleMagentization"=>irreversibleMagentization,"displacement"=> displacement)
end