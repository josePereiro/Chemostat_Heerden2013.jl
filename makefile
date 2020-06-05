# Change this to point to your env bin or alias
PYTHON3 = python3 # I used python 3.7.5
JULIA = julia  --project  # I used julia 1.1.0 (2019-01-21) (do not forget the --project)

.PHONY: clean all 
all: \
	heerden_data \
	iJR904

clear: \
	heerden_data_clear \
	iJR904_clear


###############################################
# HeerdenData
###############################################
# Prepare data taken from Heerden2013 https://doi.org/10.1186/1475-2859-12-80.
# for modeling
data/processed/heerden2013___data/heerden2013___cont_cult_data.tsv: scripts/HeerdenData/1_convert_cont_cul_data.jl
	${JULIA} $^

heerden_cont_cul_data: \
	data/processed/heerden2013___data/heerden2013___cont_cult_data.tsv



heerden_data: \
	heerden_cont_cul_data

heerden_data_clear:
	rm -fr data/processed/heerden2013___data




###############################################
# iJR904
###############################################
# model mat file
data/processed/iJR904/iJR904.mat: scripts/iJR904/0_make_mat_file.py
	$(PYTHON3) $^

iJR904_mat_file: data/processed/iJR904/iJR904.mat



data/processed/iJR904/iJR904___base_model.jls: scripts/iJR904/1_prepare_base_model.jl
	$(JULIA) $^

iJR904_base_model: data/processed/iJR904/iJR904___base_model.jls



iJR904: \
	iJR904_mat_file \
	iJR904_base_model

iJR904_clear: 
	rm -fr data/processed/iJR904