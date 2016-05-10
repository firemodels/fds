
load  = {}
unload = {}

function load.slice(matchFunc)
    local nslices = 0
    for key,value in pairs(sliceinfo) do
        if (matchFunc(value)) then
            load.datafile(value.file)
            nslices = nslices + 1
        end
    end
    if (nslices > 0) then
        print(nslices .. " slice loaded")
        return 0
    end
    error("No matching slices were found.")
end

function load.vslice(matchFunc)
    for key,value in pairs(sliceinfo) do
        if (matchFunc(value)) then
            return load.vdatafile(value.file)
        end
    end
end

function load.namedslice(name)
    return load.slice(function(slice)
        return (slice.label == name)
    end)
end

function load.namedvslice(name)
    return load.vslice(function(slice)
        return (slice.label == name)
    end)
end

function load.slice_std(slice_type, axis, distance)
    return loadslice_std(slice_type, axis, distance)
end

function load.datafile(filename)
    local errorcode = loaddatafile(filename)
    if errorcode == 1 then
        error("load.datafile: could not load " .. filename)
    end
    assert(errorcode == 0, string.format("loaddatafile errorcode: %d\n",errorcode))
    return errorcode
end

function load.vdatafile(filename)
    local errorcode = loadvdatafile(filename)
    if errorcode == 1 then
        error("load.vdatafile: could not load " .. filename)
    end
    assert(errorcode == 0, string.format("loadvdatafile errorcode: %d\n",errorcode))
    return errorcode
end

function load.tour(name)
    local errorcode = loadtour(name)
    if errorcode == 1 then
        error("load.tour: could not load " .. name)
    end
    assert(errorcode == 0, string.format("loadtour errorcode: %d\n",errorcode))
    return errorcode
end

function unload.all()
    unloadall()
end

function unload.tour()
    unloadtour()
end
