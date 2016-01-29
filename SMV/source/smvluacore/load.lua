
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
end

function load.vdatafile(filename)
    local errorcode = loadvdatafile(filename)
    assert(errorcode == 0, string.format("loadvdatafile errorcode: %d\n",errorcode))
end

function unload.all()
    unloadall()
end
