## -------------------------------------------------------------------
function load_model(modelkey::String; uncompress = true)
    file = procdir("base_models.bson")
    models = ChU.load_data(file; verbose = false);
    model_dict = models[modelkey]
    model = ChU.MetNet(;model_dict...) 
    return uncompress ? ChU.uncompressed_model(model) : model
end

function load_model(modelkey::String, exp::Int; uncompress = true)
    file = procdir("base_models.bson")
    models = ChU.load_data(file; verbose = false);
    model_dict = models[modelkey][exp]
    model = ChU.MetNet(;model_dict...) 
    return uncompress ? ChU.uncompressed_model(model) : model
end