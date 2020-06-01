# Change this to point to your env bin or alias
PYTHON3 = python3 # I used python 3.7.5
JULIA = julia  --project  # I used julia 1.1.0 (2019-01-21) (do not forget the --project)

.PHONY: clean all 
all: \
	iJR904


###############################################
# iJR904
###############################################
# model mat file
data/processed/iJR904/iJR904.mat: scripts/iJR904/0_make_mat_file.py
	$(PYTHON3) $^

iJR904_mat_file: data/processed/iJR904/iJR904.mat

iJR904: \
	iJR904_mat_file