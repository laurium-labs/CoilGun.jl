include("default_components.jl")
include("results.jl")
function dictionary_api()
    extract_results(solveScenario(default_scenario))
end