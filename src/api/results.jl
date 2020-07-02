using CoilGun: Scenario
function extract_results(sln, scenario::Scenario)
    velocity = sln[3,:]
    magnetization= sln[1,:]
    irreversibleMagentization=sln[4,:]
    displacement = sln[2,:]

    result_dict = Dict("endTime" =>scenario.endTime,"velocity" => velocity, "magnetization" => magnetization, "irreversibleMagentization"=>irreversibleMagentization,"displacement"=> displacement)
end