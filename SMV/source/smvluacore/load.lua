function loadnamedslice(name)
    for key,value in pairs(sliceinfo) do
        if (value.label == name) then
            loaddatafile(value.file)
        end
    end
end
