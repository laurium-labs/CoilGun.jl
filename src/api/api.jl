include("default_components.jl")
include("results.jl")
using JSON
function apply_json_mutations(object::Any, json::Dict)
    # This function recreates `object` with any mutations from `json`
    new_fields = map(fieldnames(typeof(object))) do field
        orig_field_value = getfield(object, field)
        if haskey(json, string(field))
            if orig_field_value isa Number
                return typeof(orig_field_value)(json[string(field)])
            else
                return apply_json_mutations(orig_field_value, json[string(field)])
            end
        else
            return orig_field_value
        end
    end
    return typeof(object)(new_fields...)
end


function get_json_description(object::Any)
    json = Dict()
    fields = fieldnames(typeof(object))
    for field in fields
        orig_field_value = getfield(object, field)
        if typeof(orig_field_value) <: Unitful.AbstractQuantity
            string_rep = string(orig_field_value)
            unit_string = string_rep[findfirst(isequal(' '), string_rep) + 1: end]
            json[string(field)] = Dict("val" => orig_field_value.val, "unit" => unit_string)
        elseif typeof(orig_field_value) <:Int
            string_rep = string(orig_field_value)
            json[string(field)] = Dict("val" => orig_field_value, "unit"=>"")
        elseif !(orig_field_value isa Function)
            json[string(field)] = get_json_description(orig_field_value)
        end
    end
    return json
end

function dictionary_api(dictionary_Coil_specification::Dict)::Dict
     ui_scenario=apply_json_mutations(CoilGunDefaults.default_scenario, dictionary_Coil_specification)
     scenario = CoilGunDefaults.transform_scenario(ui_scenario)
    #  @show scenario
     extract_results(CoilGunDefaults.solve_scenario(scenario), scenario)
end
function get_default_scenario_json()::String
    return JSON.json(get_json_description(CoilGunDefaults.default_scenario))
end