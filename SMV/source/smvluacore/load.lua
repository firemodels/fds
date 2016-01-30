
load  = {}
unload = {}
function loadnamedslice(name)
    for key,value in pairs(sliceinfo) do
        if (value.label == name) then
            loaddatafile(value.file)
        end
    end
end

function load.datafile(filename)
    local errorcode = loaddatafile(filename)
    assert(errorcode == 0, string.format("loaddatafile errorcode: %d\n",errorcode))
    return errorcode
end

function load.vdatafile(filename)
    local errorcode = loadvdatafile(filename)
    assert(errorcode == 0, string.format("loadvdatafile errorcode: %d\n",errorcode))
    return errorcode
end

function load.tour(name)
    local errorcode = loadtour(name)
    assert(errorcode == 0, string.format("loadtour errorcode: %d\n",errorcode))
    return errorcode
end

function unload.all()
    unloadall()
end

function unload.tour()
    unloadtour()
end
