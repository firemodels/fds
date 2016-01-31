--- @module smv
local smv = {}
require "bounds"
require "clipping"
require "load"
require "render"
require "view"
require "tour"

-- set the defaults for renderall
render_startframe = 0
render_skipframe = 1

smv.getfinalframe = get_nglobal_times-1
getfinalframe = smv.getfinalframe

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
model = {}
_model = {
    chid = {
        get = function()
            return fdsprefix
        end,
        set = function()
            error("model.chid is read-only")
        end
    },
    slices = {
        get = function()
            return sliceinfo
        end,
        set = function()
            error("model.slices is read-only")
        end,
        -- __len = function()
        --     error("len called slices")
        -- end
    },
    -- TODO: provide this by overriding the len operator
    nslices = {
        get = function()
            return #model.slices + 1
        end,
        set = function()
            error("model.nslices is read-only")
        end,
    },
    meshes = {
        get = function()
            -- this relies on initsmvdata being called first
            return meshinfo
        end,
        set = function()
            error("model.meshes is read-only")
        end
    },
    nmeshes = {
        get = function()
            return #model.meshes + 1
        end,
        set = function()
            error("model.nmeshes is read-only")
        end,
    }
}
local model_mt = {
   -- get method
   __index = function (t,k)
       if type(_model[k]) == "function" then
           return _model[k]
       else
           return _model[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _model[k].set(v)
   end
}
setmetatable(model, model_mt)

timebar = {}
_timebar = {
    visibility = {
        get = function()
            return gettimebarvisibility()
        end,
        set = function(v)
            return settimebarvisibility(v)
        end,
        -- toggle = function ()
        --     timebar.visibility = not timebar.visibility
        -- end
    },
}
local timebar_mt = {
   -- get method
   __index = function (t,k)
       if type(_timebar[k]) == "function" then
           return _timebar[k]
       else
           return _timebar[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _timebar[k].set(v)
   end
}
setmetatable(timebar, timebar_mt)

return smv
