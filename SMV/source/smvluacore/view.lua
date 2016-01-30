
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
            get = function ()
                return getcolorbarindex()
            end,
            set = function (v)
                setcolorbarindex(v)
            end
        }
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
       assert(_view[k], "_view." .. tostring(k) .. " does not exist.")
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

window = {}
window.size = function(width, height)
    setwindowsize(width, height)
end
