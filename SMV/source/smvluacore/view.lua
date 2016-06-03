


view = {colorbar = {}, blockages = {}, color = {}}
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
    color2bar = {
        get = function ()
            return get_color2bar_colors()
        end,
        set = function (colors)
            print("setting color2bar colors")
            return set_color2bar_colors(#colors, colors)
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
    texture_flag = {
        get = function ()
            return nil -- getcolorbarflip()
        end,
        set = function (v)
            return nil -- setcolorbarflip(v)
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
    },
    colors = {
        get = function()
            return get_colorbar_colors()
        end,
        set = function(colors)
            print("setting colorbar colors")
            return set_colorbar_colors(#colors, colors)
        end
    }
}
setmetatable(view.colorbar, colorbar_mt)

local color_mt = {
   -- get method
   __index = function (t,k)
       if type(_color[k]) == "function" then
           return _color[k]
       else
           return _color[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _color[k].set(v)
   end
}
_color = {
    ambientlight = {
        -- get = function ()
        --     return getcolorbarflip()
        -- end,
        set = function (v)
            return set_ambientlight(v.r, v.g, v.b)
        end
    },
    backgroundcolor = {
        -- get = function ()
        --     return getcolorbarindex()
        -- end,
        set = function (v)
            set_backgroundcolor(v.r, v.g, v.b)
        end
    },
    blockcolor = {
        -- get = function()
        --     return get_colorbar_visibility()
        -- end,
        set = function(v)
            return set_blockcolor(v.r, v.g, v.b)
        end
    },
    blockshininess = {
        -- get = function ()
        --     return getcolorbarindex()
        -- end,
        set = function (v)
            return set_blockshininess(v)
        end
    },
    blockspecular = {
        -- get = function()
        --     return get_colorbar_visibility()
        -- end,
        set = function(v)
            return set_blockspecular(v.r, v.g, v.b)
        end
    },
    boundcolor = {
        -- get = function()
        --     return get_colorbar_visibility()
        -- end,
        set = function(v)
            return set_boundcolor(v.r, v.g, v.b)
        end
    },
    diffuselight = {
        -- get = function()
        --     return get_colorbar_visibility()
        -- end,
        set = function(v)
            return set_diffuselight(v.r, v.g, v.b)
        end
    }
}
setmetatable(view.color, color_mt)

local blockages_mt = {
   -- get method
   __index = function (t,k)
       if type(_blockages[k]) == "function" then
           return _blockages[k]
       else
           return _blockages[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _blockages[k].set(v)
   end
}

-- View Method
-- 1 - Defined in input file
-- 2 - Solid
-- 3 - Outine only
-- 4 - Outline added
-- 5 - Hidden
local viewMethodTable = {
    as_input = 1,
    solid = 2,
    outline_only = 3,
    outline_added = 4,
    hidden = 5
}
local outlineColorTable = {
    blockage = 1,
    foreground = 2
}
local locationViewTable = {
    grid = 1,
    exact = 2,
    cad = 3
}
local convertTo = function(table, strValue)
    local v = table[strValue]
    if v == nil then error("invalid view method: " .. strValue)
        else return v
        end
    end
local convertFrom = function(table, intValue)
    for key,value in pairs(table) do
        if (value == intValue) then
            return key
        end
    end
    error("invalid view method value: " .. intValue)
    end

_blockages = {
    method = {
        -- get
        -- blockage_view_method(int setting)
        set = function(setting)
            if (type(setting) == "string") then
                blockage_view_method(convertTo(viewMethodTable, setting))
            else
                error("view.blockages.method expected string but got " .. type(setting))
            end
            end
    },
    outline_color = {
        -- get
        -- blockage_view_method(int setting)
        set = function(setting)
            if (type(setting) == "string") then
                blockage_outline_color(convertTo(outlineColorTable, setting))
            else
                error("view.blockages.outline_color expected string but got " .. type(setting))
            end
            end
    },
    locations = {
        -- get
        -- blockage_view_method(int setting)
        set = function(setting)
            if (type(setting) == "string") then
                blockage_locations(convertTo(locationViewTable, setting))
            else
                error("view.blockages.locations expected string but got " .. type(setting))
            end
            end
    }
}
setmetatable(view.blockages, blockages_mt)

window = {}
window.size = function(width, height)
    setwindowsize(width, height)
end
