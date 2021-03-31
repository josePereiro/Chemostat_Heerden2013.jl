## -------------------------------------------------------------------
function load_model(modelkey::String; uncompress = true)
    BASE_MODELS = ChU.load_data(BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS[modelkey]
    model = ChU.MetNet(;model_dict...) 
    return uncompress ? ChU.uncompressed_model(model) : model
end

function load_model(modelkey::String, exp::Int; uncompress = true)
    BASE_MODELS = ChU.load_data(BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS[modelkey][exp]
    model = ChU.MetNet(;model_dict...) 
    return uncompress ? ChU.uncompressed_model(model) : model
end