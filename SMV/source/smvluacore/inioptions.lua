
local inioptions = {}

function length(table)
    local i = 0
    for key, val in pairs(table) do i = i + 1 end
    return i
end

function mockExec(name)
    print("Executing: ", name)
end

options =
    { AMBIENTLIGHT = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting ambient light to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.ambientlight = {r = r, g = g, b = b} end
        end
    , BACKGROUNDCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting background color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.backgroundcolor = {r = r, g = g, b = b} end
        end
    , BLOCKCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting block color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.blockcolor = {r = r, g = g, b = b} end
        end
    , BLOCKSHININESS = function(opt)
        local v = opt.argLines[1][1]
        print("Setting block shininess to " .. v)
        return function () view.color.blockshininess = v end
        end
    , BLOCKSPECULAR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting block specular to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.blockspecular = {r = r, g = g, b = b} end
        end
    , BOUNDCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting bound color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.boundcolor = {r = r, g = g, b = b} end
        end
    , COLORBAR = function(opt)
        local ncolors = opt.argLines[1][1]
        local texture_flag = opt.argLines[1][2]
        local contour_value = opt.argLines[1][3]
        local r, g, b
        local colors = {}
        for i=2,ncolors+1,1 do
            r = opt.argLines[i][1]
            g = opt.argLines[i][2]
            b = opt.argLines[i][3]
            table.insert(colors,{r=r,g=g,b=b})
        end
        return function ()
            print()
            -- view.colorbar.colors = colors
            -- view.colorbar.texture_flag = (textureflag ~= 0)
            -- view.colorbar.index = contour_value
            end
        end
    , COLOR2BAR = function(opt)
        local ncolors = opt.argLines[1][1]
        local r, g, b
        local colors = {}
        for i=2,ncolors+1,1 do
            r = opt.argLines[i][1]
            g = opt.argLines[i][2]
            b = opt.argLines[i][3]
            table.insert(colors,{r=r,g=g,b=b})
        end
        return function ()
            -- view.color2bar.colors = colors
            print()
            end
        end
    , COLORBAR_FLIP = function(opt)
        local setting = opt.argLines[1][1]
        return function ()
            if (setting == 1) then
                setcolorbarflip(true)
            elseif (setting == 0) then
                setcolorbarflip(false)
            else
                error("invalid colorbar setting")
            end
            end
        end
    , DIFFUSELIGHT = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting bound color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.diffuselight = {r = r, g = g, b = b} end
        end
    , DIRECTIONCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_directioncolor(r,g,b) end
        end
    , FLIP = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_flip(setting) end
        end
    , FOREGROUNDCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_foregroundcolor(r,g,b) end
        end
    , HEATOFFCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_heatoffcolor(r,g,b) end
        end
    , HEATONCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_heatoncolor(r,g,b) end
        end
    , ISOCOLORS = function(opt)
        local shininess = opt.argLines[1][1]
        local default_opaqueness = opt.argLines[1][2]
        local specular = {}
        specular.r = opt.argLines[2][1]
        specular.g = opt.argLines[2][2]
        specular.b = opt.argLines[2][3]
        local nlevels = opt.argLines[3][1]
        local i
        local colors = {}
        for i=1,nlevels,1 do
            colors[i] = {r=opt.argLines[3+i][1],
                         g=opt.argLines[3+i][2],
                         b=opt.argLines[3+i][3]}
        end
        return function ()
            set_isocolors(shininess, default_opaqueness,
                          specular, colors)
        end
        end
    , COLORTABLE = function(opt)
        local ncolors = opt.argLines[1][1]
        local i
        local colors = {}
        for i=1,ncolors,1 do
            print(opt.argLines[i+1][1],opt.argLines[i+1][2]
                ,opt.argLines[i+1][3]
                ,opt.argLines[i+1][4]
                , opt.argLines[i+1][5]
                , opt.argLines[i+1][6]
                )
            colors[opt.argLines[i+1][6]]
                 = {r=opt.argLines[1+i][1],
                         g=opt.argLines[1+i][2],
                         b=opt.argLines[1+i][3]}
        end
        return function ()
            set_colortable(colors)
        end
        end
    , LIGHT0 = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_light0(setting) end
        end
    , LIGHT1 = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_light1(setting) end
        end
    , LIGHTMODELLOCALVIEWER = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_light1(setting) end
        end
    , LIGHTMODELSEPARATESPECULARCOLOR = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_lightmodelseparatespecularcolor(setting) end
        end
    , LIGHTPOS0 = function(opt)
        local x = opt.argLines[1][1]
        local y = opt.argLines[1][2]
        local z = opt.argLines[1][3]
        local w = opt.argLines[1][3]
        return function () set_lightpos0(x,y,z,w) end
        end
    , LIGHTPOS1 = function(opt)
        local x = opt.argLines[1][1]
        local y = opt.argLines[1][2]
        local z = opt.argLines[1][3]
        local w = opt.argLines[1][3]
        return function () set_lightpos1(x,y,z,w) end
        end
    , SENSORCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sensorcolor(r,g,b) end
        end
    , SENSORNORMCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sensornormcolor(r,g,b) end
        end
    , SETBW = function(opt)
        local geo_setting = opt.argLines[1][1]
        local data_setting = opt.argLines[1][2]
        return function () set_bw(geo_setting, data_setting) end
        end
    , SPRINKOFFCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sprinkleroffcolor(r,g,b) end
        end
    , SPRINKONCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sprinkleroncolor(r,g,b) end
        end
    , STATICPARTCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_staticpartcolor(r,g,b) end
        end
    , TIMEBARCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_timebarcolor(r,g,b) end
        end
    , VENTCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_ventcolor(r,g,b) end
        end
    , GRIDLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_gridlinewidth(v) end
        end
    , ISOLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_isolinewidth(v) end
        end
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
                assert(testArgType(argType, arg), "Argument " .. i .. " is \""
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

function inioptions.validate(v, execute)
    local matchedOption = options[v.name]
    if (matchedOption)
        then
            local execFunc = matchedOption(v)
            if execute then
                return execFunc()
            else
                return execFun
            end
        else
            error("Option " .. v.name .. " not recognised.")
    end
end

return inioptions
