--- @module smv
local smv = {}
require "clipping"
-- set the defaults for renderall
render_startframe = 0
render_skipframe = 1

function smv.render(...)
    local fstArg = select(1,...)
    if type(fstArg) == "string"
        then
            if select("#",...) == 1
            -- arg is simply a string to use as the basename
            then renderC(fstArg) -- TODO: this is unncessary, as we could
                                 -- just send all the arguments to string.format
                                 -- and it will still be happy
            -- arg is a string to be interpolated and a set of variables
            -- TODO: this doesn't work without being wrapped in a function
            -- because the arguments are strictly evaluated. This is really only
            -- relant when using renderall.
            else renderC(string.format(...))
            end
    -- arg is a function that returns a string to use as the basename
    elseif type(fstArg) == "function"
        then renderC(fstArg()) -- by only using the fstArg we strip any
                               -- extraneous arguments
    -- there is no argument and we stick with the default naming
    elseif type(fstArg) == "nil"
        then renderC(...)
    else error(type(fstArg) .. " is an invalid argument to render()")
    end
    -- TODO: return information about the file produced or something
end
render = smv.render

function smv.renderstart(startframe, skipframe)
    render_startframe = startframe
    render_skipframe = skipframe
end
renderstart = smv.renderstart

smv.getfinalframe = get_nglobal_times
getfinalframe = smv.getfinalframe


-- TODO: we may need to shift this down to the c_api in order to get the best
-- performance.
function smv.rendermany(start, final, interval, ...)
    print("luascript: Rendering every " .. interval ..
        " frame(s) starting at frame " .. start .. " until " .. final .. "\n")
    local nframes = get_nglobal_times()
    local adjFinal = nframes - 1
    if final > (nframes-1) then final = final end
    for i=start,final,interval do
        setframe(i)
        render(...)
    end
end

bounds = {}
bounds.slices = {}
-- bounds.slices["VIS_C0.9H0.1"] = {}
-- bounds.slices["VIS_C0.9H0.1"].x = {}
-- function bounds.slices["VIS_C0.9H0.1"].x.set(xMin, xMax)
function bounds.slices.set(name, min, max)
    if (min == nil)
        then
            set_slice_bound_min(name, false, 0)
        else
            set_slice_bound_min(name, true, min)
    end
    if (max == nil)
        then
            set_slice_bound_max(name, false, 0)
        else
            set_slice_bound_max(name, true, max)
    end
end

function smv.renderall(...)
    smv.rendermany(0,smv.getfinalframe(),1,...)
end
renderall = smv.renderall

function smv.loadnamedslice(name)
    for key,value in pairs(sliceinfo) do
        if (value.label == name) then
            loaddatafile(value.file)
        end
    end
end
loadnamedslice = smv.loadnamedslice

view = {colorbar = {}}
_view = {
    framenumber = {
        get =  function ()
            return getframe()
        end,
        set = function (v)
            return setframe(v)
        end
    },
    colorbar = {
        flip = {
            get = function ()
                return getcolorbarflip()
            end,
            set = function (v)
                return setcolorbarflip(v)
            end
        },
        index = {
            get = nil,
            set = function (v)
                setcolorbarindex(v)
            end
        }
    }
}
local view_mt = {
   -- get method
   __index = function (t,k)
       if type(_view[k]) == "function" then
           return _view[k]
       else
           return _view[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _view[k].set(v)
   end
}
setmetatable(view, view_mt)

local colorbar_mt = {
   -- get method
   __index = function (t,k)
       if type(_view.colorbar[k]) == "function" then
           return _view.colorbar[k]
       else
           return _view.colorbar[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _view.colorbar[k].set(v)
   end
}
setmetatable(view.colorbar, colorbar_mt)

function smv.settimeend()
    nframes = get_nglobal_times()
    setframe(nframes-1)
end
settimeend = smv.settimeend

function smv.togglecolorbarflip()
    setcolorbarflip(1-getcolorbarflip())
end
togglecolorbarflip = smv.togglecolorbarflip

function smv.colorbarnormal()
    setcolorbarflip(1)
end
colorbarnormal = smv.colorbarnormal

return smv
