


clip = {}

_clip = {
    mode = {
        get = function ()
            return getclippingmode()
        end,
        set = function (v)
            return set_clipping_mode(v)
        end
    },
    set = function(xMin, xMax, yMin, yMax, zMin, zMax)
        clip.x.set(xMin, xMax)
        clip.y.set(yMin, yMax)
        clip.z.set(zMin, zMax)
    end
}

-- the real table value
clip._x = {
   set = function (min, max)
       clip.x.max = min
       clip.x.max = max
   end,
   min = {
       set = function (v)
           if v ~= nil
               then set_sceneclip_x_min(true, v)
               else set_sceneclip_x_min(false, 0)
           end
       end
   },
   max = {
       set = function (v)
           if v ~= nil
               then set_sceneclip_x_max(true, v)
               else set_sceneclip_x_max(false, 0)
           end
       end,
   }

}
clip.x = {} -- the proxy
local x_mt = {
   -- get method
   __index = function (t,k)
       if type(clip._x[k]) == "function" then
           return clip._x[k]
       else
           return clip._x[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       clip._x[k].set(v)
   end
}
setmetatable(clip.x, x_mt)

-- the real table value
clip._y = {
   set = function (min, max)
       clip.y.max = min
       clip.y.max = max
   end,
   min = {
       set = function (v)
           if v ~= nil
               then set_sceneclip_y_min(true, v)
               else set_sceneclip_y_min(false, 0)
           end
       end
   },
   max = {
       set = function (v)
           if v ~= nil
               then set_sceneclip_y_max(true, v)
               else set_sceneclip_y_max(false, 0)
           end
       end,
   }

}
clip.y = {} -- the proxy
local y_mt = {
   -- get method
   __index = function (t,k)
       if type(clip._y[k]) == "function" then
           return clip._y[k]
       else
           return clip._y[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       clip._y[k].set(v)
   end
}
setmetatable(clip.y, y_mt)

-- the real table value
clip._z = {
   set = function (min, max)
       clip.z.max = min
       clip.z.max = max
   end,
   min = {
       set = function (v)
           if v ~= nil
               then set_sceneclip_z_min(true, v)
               else set_sceneclip_z_min(false, 0)
           end
       end
   },
   max = {
       set = function (v)
           if v ~= nil
               then set_sceneclip_z_max(true, v)
               else set_sceneclip_z_max(false, 0)
           end
       end,
   }

}
clip.z = {} -- the proxy
local z_mt = {
   -- get method
   __index = function (t,k)
       if type(clip._z[k]) == "function" then
           return clip._z[k]
       else
           return clip._z[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       clip._z[k].set(v)
   end
}
setmetatable(clip.z, z_mt)

local clip_mt = {
    -- get method
    __index = function (t,k)
        if type(_clip[k]) == "function" then
            return _clip[k]
        else
            return _clip[k].get()
        end
    end,
    -- set method
    __newindex = function (t,k,v)
        print("_clip", _clip)
        print("k", k)
        print("_clip[k]", _clip[k])
        return _clip[k].set(v)
    end
}
setmetatable(clip, clip_mt)
