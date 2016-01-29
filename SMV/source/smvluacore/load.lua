
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
    loaddatafile(filename)
end

function load.vdatafile(filename)
    loadvdatafile(filename)
end

function unload.all()
    unloadall()
end
