include("default_components.jl")
include("results.jl")
# function apply_json_mutations(object::Any, json::Dict)
#     # This function recreates `object` with any mutations from `json`
#      new_fields = map(fieldnames(typeof(object))) do field
#         orig_field_value = getfield(object, field)
#         if haskey(json, string(field))
#             if orig_field_value isa Number
#                 return typeof(orig_field_value)(json[string(field)])
#             else
#                 return apply_json_mutations(orig_field_value, json[string(field)])
#             end
#         else
#             return orig_field_value
#         end
#     end
#     return typeof(object)(new_fields...)
# end
function dictionary_api()#dictionary_Coil_specification::Dict)::Dict
    # scenario=apply_json_mutations(CoilGunDefaults.default_scenario)
     extract_results(CoilGunDefaults.solve_senario(CoilGunDefaults.default_scenario))
    # x=rand(40)
    # return string(x)
end