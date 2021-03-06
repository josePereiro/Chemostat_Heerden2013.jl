# The intakes bounds of the network are determined by the 
# medium concentration in the Chemostat model (see Cossios paper)
# This is a base medium for modeling
function load_base_intake_info()
    base_intake_info = Dict(
        "EX_glc_LPAREN_e_RPAREN_" => Dict("c"=>  100.0, "lb"=> -100.0),
        "EX_nh4_LPAREN_e_RPAREN_" => Dict("c"=>99999.0, "lb"=> -100.0),
        "EX_o2_LPAREN_e_RPAREN_"  => Dict("c"=>99999.0, "lb"=> -100.0),
        "EX_pi_LPAREN_e_RPAREN_"  => Dict("c"=>99999.0, "lb"=> -100.0),
        "EX_so4_LPAREN_e_RPAREN_" => Dict("c"=>99999.0, "lb"=> -100.0),
        "EX_h2o_LPAREN_e_RPAREN_" => Dict("c"=>99999.0, "lb"=> -100.0),
        "EX_h_LPAREN_e_RPAREN_" => Dict("c"=>99999.0, "lb"=> -100.0),
    )
end