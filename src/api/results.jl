using CoilGun: Scenario
function extract_results(sln, scenario::Scenario)
    velocity = sln[3,:]
    magnetization= sln[1,:]
    irreversibleMagentization=sln[4,:]
    displacement = sln[2,:]
    endTime = scenario.endTime
    result_dict = Dict( "velocity" => velocity, "magnetization" => magnetization, "irreversibleMagentization"=>irreversibleMagentization,"displacement"=> displacement,"endTime"=>endTime)
end