# Here I include some maps between the experimental
# data ids and the model ids

######################
# maps between Heerden2013 https://doi.org/10.1186/1475-2859-12-80 
# data and the model
#####################

# Mets map
Hd_mets_map = Dict();
Hd_mets_map["GLC"] = "glc_DASH_D_e";
Hd_mets_map["SA"] = "succ_e";
Hd_mets_map["AcA"] = "ac_e";
Hd_mets_map["FA"] = "for_e";
Hd_mets_map["MA"] = "mal_DASH_D_e";
for (k, v) in Hd_mets_map
    Hd_mets_map[v] = k
end

# Map for Heerden2013 https://doi.org/10.1186/1475-2859-12-80. Table 2
Hd_to_inactivate_map = Dict(
    "2-ketobutyrate formate lyase"=>"OBTFL",
    "Acetate kinase"=>"ACKr",
    "Alcohol dehydrogenase"=>"ALCD19",
    "Citrate lyase"=>"CITL",
    "Lactate dehydrogenase"=>"LDH_D",
    "Methylglyoxal synthase"=>"MGSA",
    "NAD+-linked malic enzyme"=>"ME1",
    "Phosphotransacetylase"=>"PTAr",
    "Pyruvate formate lyase"=>"PFL",
    "Pyruvate oxidase"=>"POX",
    # "Threonine decarboxylase"=> "" # not found in model TODO: Find it
)

# the model has no way to simulate overexpression
Hd_to_activate_map = Dict(
    "PEP carboxikinase"=>["PPCK"]
)



###################
# enzymatic costs
# from Beg, (2007) https://doi.org/10.1073/pnas.0609845104.
###################

# A map between model ids and the reactions reported in Beg2007
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
"DNA ligase (NAD+)"=>[""])

# base model exch met map
# A quick way to get exchages from mets and te other way around
const exch_met_map = Dict()
function _load_exch_met_map()
    !isfile(EXCH_MET_MAP_FILE) && return exch_met_map
    empty!(exch_met_map)
    merge!(exch_met_map, load_data(EXCH_MET_MAP_FILE; verbose = false))
    return exch_met_map
end
