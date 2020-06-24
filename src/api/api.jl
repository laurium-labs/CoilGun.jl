include("default_components.jl")
include("results.jl")
function dictionary_api()
    extract_results(CoilGunDefaults.solve_senario(CoilGunDefaults.default_scenario))
    # x=rand(40)
    # return string(x)
end