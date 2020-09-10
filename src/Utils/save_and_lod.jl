function load_data(src_file; verbose = true)
    data = wload(src_file)[DATA_KEY]
    verbose && println(relpath(src_file), " loaded!!!, size: ", filesize(src_file), " bytes")
    return data
end

function save_data(src_file, data; verbose = true)
    data = tagsave(src_file, Dict(DATA_KEY => data))
    verbose && println(relpath(src_file), " saved!!!, size: ", filesize(src_file), " bytes")
    return data
end
