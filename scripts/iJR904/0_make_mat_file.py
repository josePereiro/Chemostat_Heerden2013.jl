# -*- coding: utf-8 -*-
import os
import errno
import cobra

# TODO find a way to do this in COBRA.jl

# +
# Here we take iJR904 model and make a mat_file from it.
# -

# %cd /Users/Pereiro/.julia/dev/Chemostat_Heerden2013

# ### Meta

model_name = 'iJR904'


# ### Loading SMBL
# This is the fisrt script to be executed, so it checks for the dir
raw_gems_dir = f"./data/raw/{model_name}"
if not os.path.exists(raw_gems_dir):
    raise FileNotFoundError(
        errno.ENOENT, os.strerror(errno.ENOENT), raw_gems_dir)

print(f"Model {model_name} from (Reed et al., 2003) (download link: https://darwin.di.uminho.pt/models)")
# TODO search model download link
model_file = os.path.join(raw_gems_dir, f"{model_name}.xml")
print('Loading ', model_file, "...")
model = cobra.io.read_sbml_model(model_file)
print('Model loaded', (len(model.metabolites), len(model.reactions)))


# ### Saving mat

processed_gems_dir = f"./data/processed/{model_name}"
if not os.path.exists(processed_gems_dir):
    os.makedirs(processed_gems_dir)
    print(f"created {processed_gems_dir}")

print("Saving mat file ...")
mat_file = os.path.join(processed_gems_dir, f"{model_name}.mat")
cobra.io.save_matlab_model(model, mat_file, varname = 'model')
print(f"created {mat_file}")

# ### Checking mat model

print("Cheking mat model...")
mat_model = cobra.io.mat.load_matlab_model(mat_file);
print('Mat Model loaded', (len(mat_model.metabolites), len(mat_model.reactions)))

# +
# Cheking Dimention
model_size = (len(model.metabolites), len(model.reactions))
mat_model_size = (len(mat_model.metabolites), len(mat_model.reactions))

if not model_size == mat_model_size:
    raise Exception(f"sbml and mat model dimention missmatch!!!")
# -

# Cheking reaction bounds
for (sbml_rxn, mat_rxn) in zip(model.reactions, mat_model.reactions):
    if sbml_rxn.id != mat_rxn.id:
        raise Exception(f"ids dont match!!!, {sbml_rxn.id} != {mat_rxn.id}")
    if not sbml_rxn.bounds == mat_rxn.bounds:
        raise Exception(f"sbml and mat model dimention missmatch!!!")
