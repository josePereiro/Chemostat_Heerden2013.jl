# The intakes bounds of the network are determined by the 
# medium concentration in the Chemostat model (see Cossios paper)
# This is a base medium for modeling
base_intake_info = Dict(
    "EX_glc__D_e" => Dict("c"=>  15.0, "lb"=> -100.0),
    "EX_h2o_e" => Dict("c"=>  99999.0, "lb"=> -100.0),
    "EX_nh4_e" => Dict("c"=>  99999.0, "lb"=> -100.0),
    "EX_o2_e" => Dict("c"=>  99999.0, "lb"=> -100.0),
    "EX_pi_e" => Dict("c"=>  99999.0, "lb"=> -100.0),
    # "EX_acald_e" => Dict("c"=>  99999.0, "lb"=> -100.0),
)