using CoilGun: Scenario
function extract_results(sln)
    velocity = sln[3,:]
    magnetization= sln[1,:]
    irreversibleMagentization=sln[4,:]
    displacement = sln[2,:]

    result_dict = Dict("velocity" => velocity, "magnetization" => magnetization, "irreversibleMagentization"=>irreversibleMagentization,"displacement"=> displacement)
end