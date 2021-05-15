# Here I include some maps between the experimental
# data ids and the model ids

######################
# maps between Heerden2013 https://doi.org/10.1186/1475-2859-12-80 
# data and the model
#####################
# Mets map
function load_mets_map()   
    mets_map = Dict()
    mets_map["GLC"] = "glc_DASH_D_e"
    mets_map["SUCC"] = "succ_e"
    mets_map["FORM"] = "for_e"
    mets_map["THM"] = "thm_e"
    mets_map["NH4"] = "nh4_e"
    mets_map["CIT"] = "cit_e"
    mets_map["CO2"] = "co2_e"
    mets_map["O2"] = "o2_e"
    mets_map["AC"] = "ac_e"
    mets_map["PYR"] = "pyr_e"
    mets_map["LAC"] = "lac_DASH_D_e"
    mets_map["MAL"] = "mal_DASH_D_e";
    for (k, v) in mets_map
        mets_map[v] = k
    end
    return mets_map
end

function load_rxns_map() 
    rxns_map = Dict()
    rxns_map["D"] = "BiomassEcoli"
    rxns_map["AC"] = "EX_ac_LPAREN_e_RPAREN_"
    rxns_map["NH4"] = "EX_nh4_LPAREN_e_RPAREN_"
    rxns_map["GLC"] = "EX_glc_LPAREN_e_RPAREN_"
    rxns_map["THM"] = "EX_thm_LPAREN_e_RPAREN_"
    rxns_map["CO2"] = "EX_co2_LPAREN_e_RPAREN_"
    rxns_map["FORM"] = "EX_for_LPAREN_e_RPAREN_"
    rxns_map["CIT"] = "EX_cit_LPAREN_e_RPAREN_"
    rxns_map["SUCC"] = "EX_succ_LPAREN_e_RPAREN_"
    rxns_map["O2"] = "EX_o2_LPAREN_e_RPAREN_"
    rxns_map["MAL"] = "EX_mal_L_LPAREN_e_RPAREN_"
    rxns_map["LAC"] = "EX_lac_D_LPAREN_e_RPAREN_"
    rxns_map["PYR"] = "EX_pyr_LPAREN_e_RPAREN_"
    for (k, v) in rxns_map
        rxns_map[v] = k
    end
    return rxns_map
end

# Map for Heerden2013 https://doi.org/10.1186/1475-2859-12-80. Table 2
function load_Hd_to_inactivate_map()
    Hd_to_inactivate_map = Dict(

        # P42632
        # EC:2.3.1. 
        # 2-oxobutanoate + CoA = formate + propanoyl-CoA
        "2-ketobutyrate formate lyase"=>["OBTFL"],

        # P0A6A3
        # EC:2.7.2.1
        # acetate + ATP = acetyl phosphate + ADP 
        "Acetate kinase"=>["ACKr"],

        # (-1.0) glyald_c + (-1.0) h_c + (-1.0) nadh_c + (-0.00242) cost ==> (1.0) glyc_c + (1.0) nad_c
        "Alcohol dehydrogenase"=>["ALCD19"],

        # (-1.0) asp_DASH_L_c + (-1.0) cbp_c + (-4.4e-5) cost ==> (1.0) cbasp_c + (1.0) h_c + (1.0) pi_c
        "Aspartate aminotransferase" => ["ASPCT"],

        # (-1.0) cit_c + (-0.0031) cost ==> (1.0) ac_c + (1.0) oaa_c
        "Citrate lyase"=>["CITL"],

        # (-1.0) for_e + (-0.0031) cost ==> (1.0) for_c
        "Formate transporter" => ["FORt"],

        # (-1.0) lac_DASH_D_c + (-1.0) q8_c + (-0.0031) cost ==> (1.0) pyr_c + (1.0) q8h2_c
        "Lactate dehydrogenase"=> ["LDH_D2", "LDH_D"],

        # P0A731
        # (-1.0) dhap_c + (-0.0031) cost ==> (1.0) mthgxl_c + (1.0) pi_c
        "Methylglyoxal synthase"=>["MGSA"],

        # P26616
        # (-1.0) mal_DASH_L_c + (-1.0) nad_c + (-0.0031) cost ==> (1.0) co2_c + (1.0) nadh_c + (1.0) pyr_c
        "NAD+-linked malic enzyme"=>["ME1"],

        # P0A9M8
        # acetyl-CoA + phosphate = acetyl phosphate + CoA
        "Phosphotransacetylase"=>["PTAr"],

        # P42632
        # 2-oxobutanoate + CoA = formate + propanoyl-CoA
        "Pyruvate formate lyase"=>["PFL"],

        # (-1.0) h2o_c + (-1.0) pyr_c + (-1.0) q8_c + (-0.0031) cost ==> (1.0) ac_c + (1.0) co2_c + (1.0) q8h2_c
        "Pyruvate oxidase"=>["POX"],

        # "Threonine decarboxylase"=> "" # not found in model TODO: Find it
    )
end

# the model has no way to simulate overexpression
function load_Hd_to_activate_map()
    Hd_to_activate_map = Dict(
        "PEP carboxikinase"=>["PPCK"]
    )
end



###################
# enzymatic costs
# from Beg, (2007) https://doi.org/10.1073/pnas.0609845104.
###################

# A map between model ids and the reactions reported in Beg2007
function load_beg_rxns_map()
    beg_rxns_map = Dict("carbonyl reductase (NADPH)"=>["P5CR"],
        "alcohol dehydrogenase (NADP+)"=>["ALCD19"],
        "quinate/shikimate dehydrogenase"=>["SHK3Dr"],
        "malate dehydrogenase (decarboxylating)"=>["MDH","MDH2", "MDH3"],
        "3alpha-hydroxysteroid dehydrogenase (B-specific)"=>[""],
        "2-hydroxy-3-oxopropionate reductase"=>[""],
        "glucose dehydrogenase (acceptor)"=>["G6PDH2r", "UDPGD"],
        "cellobiose dehydrogenase (acceptor)"=>[""],
        "peroxidase"=>[""],
        "catechol 2,3-dioxygenase"=>[""],
        "arachidonate 8-lipoxygenase"=>[""],
        "calcidiol 1-monooxygenase"=>[""],
        "nitric-oxide synthase"=>[""],
        "phenylalanine 4-monooxygenase"=>[""],
        "tryptophan 5-monooxygenase"=>[""],
        "Carboxylate reductase"=>["P5CR"],
        "arsenate reductase (donor)"=>[""],
        "biliverdin reductase"=>[""],
        "15-oxoprostaglandin 13-oxidase"=>[""],
        "coproporphyrinogen oxidase"=>["CPPPGO","PPPGO"],
        "long-chain-acyl-CoA dehydrogenase"=>[""],
        "butyryl-CoA dehydrogenase"=>[""],
        "acyl-CoA dehydrogenase"=>[""],
        "L-amino-acid oxidase"=>["ASPO3","ASPO4","ASPO5","ASPO6"],
        "amine oxidase (flavin-containing)"=>["PYAM5PO"],
        "methylenetetrahydrofolate reductase [NAD(P)H]"=>["MTHFR2"],
        "formyltetrahydrofolate dehydrogenase"=>["MTHFD"],
        "sarcosine oxidase"=>[""],
        "nitrate reductase (NADH)"=>[""],
        "nitrite reductase (NO-forming)"=>[""],
        "nitrate reductase"=>["NO3R2"],
        "trypanothione-disulfide reductase"=>[""],
        "glutathione-disulfide reductase"=>[""],
        "thioredoxin-disulfide reductase"=>[""],
        "thiol oxidase"=>[""],
        "nitrate reductase (cytochrome)"=>[""],
        "aspartate carbamoyltransferase"=>["ASPCT"],
        "serine O-acetyltransferase"=>["SERAT"],
        "protein-glutamine gamma-glutamyltransferase"=>[""],
        "gamma-glutamyltransferase"=>["CRNBTCT"],
        "citrate (Si)-synthase"=>["CS"],
        "kaempferol 3-O-galactosyltransferase"=>[""],
        "NAD+ ADP-ribosyltransferase"=>["NNDMBRT"],
        "di-trans,poly-cis-decaprenylcistransferase"=>["ACGAMT"],
        "cystathionine gamma-synthase"=>[""],
        "adenosine kinase"=>["ADNK1"],
        "glycerate kinase"=>["GLYCK"],
        "galactokinase"=>["GALKr"],
        "[pyruvate dehydrogenase (acetyl-transferring)] kinase"=>[""],
        "guanylate kinase"=>["GK1"],
        "FMN adenylyltransferase"=>["FMNAT"],
        "tRNA adenylyltransferase"=>[""],
        "aryl sulfotransferase"=>[""],
        "aminoacyl-tRNA hydrolase"=>[""],
        "carboxymethylenebutenolidase"=>[""],
        "ubiquitin thiolesterase"=>[""],
        "fructose-bisphosphatase"=>["FBP","FBA"],
        "[phosphorylase] phosphatase"=>[""],
        "phosphoglycolate phosphatase"=>["PGLYCP"],
        "protein-tyrosine-phosphatase"=>[""],
        "inositol-polyphosphate 5-phosphatase"=>["MI1PP"],
        "3',5'-cyclic-GMP phosphodiesterase"=>[""],
        "beta-glucosidase"=>["MLTG1","MLTG2","MLTG3","MLTG4","MLTG5"],
        "beta-glucuronidase"=>[""],
        "glucosylceramidase"=>[""],
        "cyclomaltodextrinase"=>[""],
        "alpha-N-arabinofuranosidase"=>[""],
        "purine nucleosidase"=>["AMPN","AHCYSNS","CMPN","MTAN","NMNN"],
        "rRNA N-glycosylase"=>[""],
        "NAD+ nucleosidase"=>[""],
        "Xaa-Pro aminopeptidase"=>[""],
        "dipeptidyl-peptidase I"=>[""],
        "peptidyl-dipeptidase A"=>[""],
        "coagulation factor Xa"=>[""],
        "t-Plasminogen activator"=>[""],
        "cathepsin B"=>[""],
        "envelysin"=>[""],
        "amidase"=>["","GSPMDA","NMNDA","NNAM"],
        "formamidase"=>[""],
        "arginase"=>[""],
        "guanidinoacetase"=>["GUAD"],
        "apyrase"=>[""],
        "phloretin hydrolase"=>[""],
        "Orotidine-5'-phosphate decarboxylase"=>["OMPDC"],
        "4-Hydroxybenzoate decarboxylase"=>["OPHBDC"],
        "Threonine aldolase"=>["THRAr"],
        "enoyl-CoA hydratase"=>[""],
        "Uroporphyrinogen-III synthase"=>["UPP3S"],
        "dihydroxy-acid dehydratase"=>["DHAD1","DHAD2"],
        "pectin lyase"=>[""],
        "DNA-(apurinic or apyrimidinic site) lyase"=>[""],
        "lactoylglutathione lyase"=>["LGTHL"],
        "guanylate cyclase"=>[""],
        "dTDP-4-dehydrorhamnose 3,5-epimerase"=>["TDPDRE"],
        "UDP-glucose 4-epimerase"=>["UDPG4E"],
        "Triose-phosphate isomerase"=>["TPI"],
        "steroid DELTA-isomerase"=>[""],
        "dodecenoyl-CoA isomerase"=>[""],
        "Glutamate-1-semialdehyde 2,1-aminomutase"=>[""],
        "Chalcone isomerase"=>[""],
        "Chloromuconate cycloisomerase"=>[""],
        "Tyrosine-tRNA ligase"=>[""],
        "Threonine-tRNA ligase"=>[""],
        "Isoleucine-tRNA ligase"=>[""],
        "Lysine-tRNA ligase"=>[""],
        "formate-tetrahydrofolate ligase"=>[""],
        "Adenylosuccinate synthase"=>["ADSS"],
        "DNA ligase (NAD+)"=>[""]
    )
end

# base model exch met map
# A quick way to get exchages from mets and te other way around
function load_exch_met_map()
    datfile = proddir("exch_met_map.bson")
    !isfile(datfile) && return Dict()
    return load_data(datfile; verbose = false)
end

function load_krebs_iders()
    krebs_iders = ["SUCD1", "SUCOAS", "AKGDH", "ICDHyr", 
        "ACONT", "CS", "MDH", "FUM", "MALS", "ICL"
    ]
end

function load_inners_idermap()
    kreps_idermap = Dict(
        "SUCD1" => ["SUCD1i"], 
        "SUCOAS" => ["SUCOAS_fwd", "SUCOAS_bkwd"], 
        "AKGDH" => ["AKGDH"],
        "ICDHyr" => ["ICDHyr_fwd", "ICDHyr_bkwd"],
        "ACONT" => ["ACONT_bkwd", "ACONT_fwd"],
        "CS" => ["CS"],
        "MDH" => ["MDH_fwd", "MDH_bkwd"],
        "FUM" => ["FUM_fwd", "FUM_bkwd"],
        "MALS" => ["MALS"],
        "ICL" => ["ICL"]
    )
end