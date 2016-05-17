


view = {colorbar = {}}
_view = {
    -- colorbar = {
    --     get = function()
    --         return get_colorbar_visibility()
    --     end,
    --     set = function(setting)
    --         if (type(setting) == "boolean")
    --             then set_colorbar_visibility(setting)
    --         else
    --             error("the argument of set_colorbar_visibility must be a boolean, but it is a " .. type(setting))
    --         end
    --     end
    -- },
    framenumber = {
        get =  function ()
            return getframe()
        end,
        set = function (v)
            return setframe(v)
        end
    },
    viewpoint = {
        get = function ()
            return getviewpoint()
        end,
        set = function (v)
            local errorcode = setviewpoint(v)
            assert(errorcode == 0, string.format("setviewpoint errorcode: %d\n",errorcode))
            return errorcode
        end
    },
    projection_type = {
        get = function ()
            return camera_get_projection_type()
        end,
        set = function (v)
            if not (type(v) == "number" and (v == 0 or v == 1)) then
              error("projection type: " .. v .. " invalid")
            end
            local errorcode = camera_set_projection_type(v)
            assert(errorcode == 0, string.format("set_projection_type errorcode: %d\n",errorcode))
            return errorcode
        end
    }
}
local view_mt = {
   -- get method
   __index = function (t,k)
      if type(_view[k]) == "function" then
           return _view[k]
      -- elseif type(_view[k]) == "table" then
      --       return _view[k]
      else
          return _view[k].get()
      end
   end,
   -- set method
   __newindex = function (t,k,v)
       assert(_view[k], "_view." .. tostring(k) .. " does not exist.")
       _view[k].set(v)
   end
}
setmetatable(view, view_mt)

local colorbar_mt = {
   -- get method
   __index = function (t,k)
       if type(_colorbar[k]) == "function" then
           return _colorbar[k]
       else
           return _colorbar[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _colorbar[k].set(v)
   end
}
_colorbar = {
    flip = {
        get = function ()
            return getcolorbarflip()
        end,
        set = function (v)
            return setcolorbarflip(v)
        end
    },
    index = {
        get = function ()
            return getcolorbarindex()
        end,
        set = function (v)
            setcolorbarindex(v)
        end
    },
    show = {
        get = function()
            return get_colorbar_visibility()
        end,
        set = function(setting)
            if (type(setting) == "boolean")
                then set_colorbar_visibility(setting)
            else
                error("the argument of set_colorbar_visibility must be a boolean, but it is a " .. type(setting))
            end
        end
    }
}

setmetatable(view.colorbar, colorbar_mt)

window = {}
window.size = function(width, height)
    setwindowsize(width, height)
end
