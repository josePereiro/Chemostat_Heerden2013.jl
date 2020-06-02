# TODO finish this
# function add_Heerden_genes_mods(metnet::MetNet; verbose = false)

#     verbose && print("Modifying reations")
#     for (Hd_id, model_ids) in to_inactivate_map
#         # Closing formate tranportation
#         if Hd_id == "Formate transporter"
#             METS.set_L!(0.0, Hd_mets_map["FA"])
#             METS.set_V!(0.0, Hd_mets_map["FA"])
#             verbose && print(Hd_id, "blocked")
#         else
#             # Blocking inactivated reactions
#             for rxn_id in model_ids
#                 RXNS.set_pub!(0.0, rxn_id)
#                 RXNS.set_nub!(0.0, rxn_id)
#                 verbose && println(rxn_id, "blocked")
#             end
#         end
#     end

# end