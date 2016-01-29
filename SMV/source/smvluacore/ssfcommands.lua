
local ssfcommands = {}

function length(table)
    local i = 0
    for key, val in pairs(table) do i = i + 1 end
    return i
end

function mockExec(name)
    print("Executing: ", name)
end

commands =
    { CBARFLIP =    {nargs = 0, argTypes = {},
        func = function(args) return colorbarflip() end}
    , CBARNORMAL =  {nargs = 0, argTypes = {},
        func = function(args) return colorbarnormal() end}
    , EXIT =        {nargs = 0, argTypes = {},
        func = function(args) return exit() end}
    , KEYBOARD =    {nargs = 0, argTypes = {},
        func = function(args) return mockExec("KEYBOARD") end}
    , GSLICEORIEN = {nargs = 2, argTypes = {"float", "float"},
        func = function(args) return gsliceorien(args[1], args[2]) end}
    , GSLICEPOS =   {nargs = 3, argTypes = {"float", "float", "float"},
        func = function(args) return gslicepos(args[1], args[2], args[3]) end}
    , GSLICEVIEW =  {nargs = 4, argTypes = {"integer", "boolean", "boolean",
        "boolean"},
        func = function(args) return gsliceview(args[1], args[2], args[3]) end}
    , LOAD3DSMOKE = {nargs = 1, argTypes = {"string"},
        func = function(args) return load3dsmoke(args[1]) end}
    , LOADBOUNDARY = {nargs = 1, argTypes = {"string"},
        func = function(args) return loadboundary(args[1]) end}
    , LOADFILE =    {nargs = 1, argTypes = {"string"},
        func = function(args) return loaddatafile(args[1]) end}
    , LOADINIFILE = {nargs = 1, argTypes = {"string"},
        func = function(args) return loadinifile(args[1]) end}
    , LOADISO =     {nargs = 1, argTypes = {"string"},
        func = function(args) return loadiso(args[1]) end}
    , LOADPARTICLES = {nargs = 1, argTypes = {"string"},
        func = function(args) return loadparticles(args[1]) end}
    , LOADPLOT3D =  {nargs = 2, argTypes = {"int", "float"},
        func = function(args) return loadplot3d(args[1], args[2]) end}
    , LOADTOUR = {  nargs = 1, argTypes = {"string"},
        func = function(args) return loadtour(args[1]) end}
    , LOADVOLSMOKE = {nargs = 1, argTypes = {"int"},
        func = function(args) return loadvolsmoke(args[1]) end}
    , LOADVOLSMOKEFRAME = {nargs = 2, argTypes = {"int", "int"},
        func = function(args) return loadvolsmokeframe(args[1], args[2]) end}
    , LOADVFILE = {nargs = 1, argTypes = {"string"},
        func = function(args) return loadvdatafile(args[1]) end}
    , LOADVSLICE = {nargs = 3, argTypes = {"string", "int", "float"},
        func = function(args) return loadvslice(args[1], args[2], args[3]) end}
    , MAKEMOVIE = {nargs = 3, argTypes = {"string", "string", "float"},
        func = function(args) return makemovie(args[1], args[2], args[3]) end}
    , PARTCLASSCOLOR = {nargs = 1, argTypes = {"string"},
        func = function(args) return partclasscolor(args[1]) end}
    , PARTCLASSTYPE = {nargs = 1, argTypes = {"string"},
        func = function(args) return partclasstype(args[1]) end}
    , PLOT3DPROPS = {nargs = 5, argTypes = {"int", "int", "int", "int", "int",
        "float"},
        func = function(args) return plot3dprops(args[1], args[2], args[3], args[4],
            args[5]) end}
    , RENDERALL = {nargs = 0, argTypes = {},
        func = function(args) return renderall() end}
    , RENDERCLIP = {nargs = 5, argTypes = {"boolean", "int", "int", "int", "int"},
        func = function(args) return renderclip(args[1],args[2], args[3], args[4],
            args[5]) end}
    , RENDERDIR = {nargs = 1, argTypes = {{optional = true, type = "string"}},
        func = function(args) return setrenderdir(args[1]) end}
    , RENDERTYPE = {nargs = 1, argTypes = {"string"},
        func = function(args) return rendertype(args[1]) end}
    , MOVIETYPE = {nargs = 1, argTypes = {"string"},
        func = function(args) return movietype(args[1]) end}
    , RENDERSIZE = {nargs = 0, argTypes = {},
        func = function(args) return error("RENDERSIZE not implemented") end}
    , RENDERDOUBLEONCE = {nargs = 0, argTypes = {},
        func = function(args) return error("RENDERDOUBLEONCE not implemented") end}
    , RENDERONCE = {nargs = 0, argTypes = {},
        func = function(args) return render() end}
    , RENDERSTART = {nargs = 0, argTypes = {},
        func = function(args) return mockExec("command") end}
    , SCENECLIP = {nargs = 1, argTypes = {"int"},
        func = function(args) return set_sceneclip(args[1]) end}
    , SETTOURKEYFRAME = {nargs = 1, argTypes = {"float"},
        func = function(args) return settourkeyframe(args[1]) end}
    , SETTOURVIEW = {nargs = 4, argTypes = {"int", "int", "boolean", "float"},
        func = function(args) return settourview(args[1], args[2], args[3], args[4])
            end}
    , SETTIMEVAL = {nargs = 1, argTypes = {"float"},
        func = function(args) return settime(args[1]) end}
    , SETVIEWPOINT = {nargs = 1, argTypes = {"string"},
        func = function(args) return setviewpoint(args[1]) end}
    , SHOWPLOT3DDATA = {nargs = 0, argTypes = {},
        func = function(args) return error("SHOWPLOT3DDATA not implemented") end}
    , UNLOADALL = {nargs = 0, argTypes = {},
        func = function(args) return unloadall() end}
    , UNLOADTOUR = {nargs = 0, argTypes = {},
        func = function(args) return unloadtour() end}
    , VOLSMOKERENDERALL = {nargs = 0, argTypes = {},
        func = function(args) return error("VOLSMOKERENDERALL not implemented") end}
    , ISORENDERALL = {nargs = 0, argTypes = {},
        func = function(args) return error("ISORENDERALL not implemented") end}
    , XSCENECLIP = {nargs = 4, argTypes = {"boolean", "float", "boolean",
        "float"},
        func = function(args) return set_sceneclip_x(args[1], args[2], args[3],
            args[4]) end}
    , YSCENECLIP = {nargs = 4, argTypes = {"boolean", "float", "boolean",
        "float"},
        func = function(args) return set_sceneclip_y(args[1], args[2], args[3],
            args[4]) end}
    , ZSCENECLIP = {nargs = 4, argTypes = {"boolean", "float", "boolean",
        "float"},
        func = function(args) return set_sceneclip_z(args[1], args[2], args[3],
            args[4]) end}
    -- the below are not in the SMV source code, but exist in example files
    , LOADSLICEM = {nargs = 0, argTypes = {},
        func = function(args) return warning("LOADSLICEM does not exist in smokeview") end}
    }

-- an ssfType is either a string or one of the three number types:
-- boolean, int, and float. A value of 0 1 one will be classed as a boolean by
-- this function, even though it may function as an int or a float and so on.
-- There is the patter boolean > int > float, so if a float is required it is
-- perfectly fine for a boolean to be given, as boolean is simple the most
-- restrictive type it satisfies.
function isInt(n) return (n % 1) == 0 end
function isBool(n) return n == 1 or n == 0 end
function ssfNumberType(n)
    if isBool(n) then return "boolean"
    elseif isInt(n) then return "integer"
    else return "float"
    end
end
function ssfType(arg)
    local t = type(arg)
    if t == "string" then return t
    elseif t == "number"
        then return ssfNumberType(arg)
        else error("invalid type: " .. t)
    end
end

function testArgType(specArgType, arg)
    local usedArgType = ssfType(arg)
    if specArgType == "string" then return usedArgType == "string"
    elseif specArgType == "float" then return usedArgType == "boolean"
        or usedArgType == "integer" or usedArgType == "float"
    elseif specArgType == "integer" then return usedArgType == "boolean"
        or usedArgType == "integer"
    else return usedArgType == "boolean"
    end
end

-- does not fault on too many arguments being given, they are simply ignored
function testArgTypes(matchedCommand, v)
    local argTypes = matchedCommand.argTypes
    local args = v.args
    for i, argType in ipairs(argTypes) do
        local arg = args[i]
        local optional
        if type(argType) == "table"
            then optional = argType.optional
                 argType = argType.type
            else optional = false
        end
        if arg
            then
                assert(testArgType(argType, arg), "Argument 3 is \""
                .. arg .. "\" but should be of the type " .. argType)
            else
                if not optional then
                    local expected = #argTypes
                    local given = #args
                    error(v.command .. ": Insufficient arguments."
                        .. " Expected " .. expected .. " but was given "
                        .. given)
                end
        end
    end
    return true
end

function ssfcommands.validate(v, execute)
    local matchedCommand = commands[v.command]
    if (matchedCommand)
        then
            testArgTypes(matchedCommand, v)
            matchedCommand.args = v.args
            matchedCommand.command = v.command
            if execute then matchedCommand.func(v.args) end
            return matchedCommand
        else
            error("Command " .. v.command .. " not recognised.")
    end
end

return ssfcommands
