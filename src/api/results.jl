using CoilGun: Scenario
using Unitful: ustrip
function extract_results(sln, scenario::Scenario)
    velocity = sln[3,:]
    magnetization= sln[1,:]
    irreversibleMagentization=sln[4,:]
    displacement = sln[2,:]
    endTime = ustrip(scenario.endTime)
    result_dict = Dict( "time"=> sln.t, "velocity" => velocity, "magnetization" => magnetization, "irreversibleMagentization"=>irreversibleMagentization,"displacement"=> displacement)
end