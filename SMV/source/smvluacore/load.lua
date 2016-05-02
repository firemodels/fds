
load  = {}
unload = {}
function load.namedslice(name)
    for key,value in pairs(sliceinfo) do
        if (value.label == name) then
            loaddatafile(value.file)
        end
    end
end

function load.slice(...)
    return loadslice(...)
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
