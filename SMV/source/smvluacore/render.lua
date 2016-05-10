render = {movie = {}}
function _renderF(...)
    local errorcode
    local fstArg = select(1,...)
    if (...) == nil then errorcode = renderC(...)
    elseif type(fstArg) == "string" then
        if select("#",...) == 1
        -- arg is simply a string to use as the basename
        then errorcode = renderC(fstArg) -- TODO: this is unncessary, as we could
                             -- just send all the arguments to string.format
                             -- and it will still be happy
        -- arg is a string to be interpolated and a set of variables
        -- TODO: this doesn't work without being wrapped in a function
        -- because the arguments are strictly evaluated. This is really only
        -- relant when using renderall.
        else errorcode = renderC(string.format(...))
        end
    -- arg is a function that returns a string to use as the basename
    elseif type(fstArg) == "function" then errorcode = renderC(fstArg())
     -- by only using the fstArg we strip any
                               -- extraneous arguments
    -- there is no argument and we stick with the default naming
    elseif type(fstArg) == "nil" then errorcode = renderC(...)
    else error(type(fstArg) .. " is an invalid argument to render()")
    end
    if (errorcode == 0)
      then return errorcode
    else
      error("render returned errorcode: " .. errorcode)
    end
    -- TODO: return information about the file produced or something
end

_render = {
    type = {
        get = function ()
            return get_rendertype()
        end,
        set = function(v)
            assert(v == "JPG" or v == "PNG", v .. " is not a valid render type."
                   .. " Please choose from JPG or PNG")
            return set_rendertype(v)
        end
    },
    movie = {
        type = {
            get = function ()
                return get_movietype()
            end,
            set = function (v)
                assert(v == "WMV" or v == "MP4" or v == "AVI",
                       v .. " is not a valid render type. Please choose from "
                       .. "WMV, MP4, or AVI")
                return set_movietype(v)
            end
        }
    },
    dir = {
        get = function ()
            return getrenderdir();
        end,
        set = function(v)
            return setrenderdir(v)
        end
    }
}
local render_mt = {
    __call = function (t,...)
        return _renderF(...)
    end,
   -- get method
   __index = function (t,k)
       if type(_render[k]) == "function" then
           return _render[k]
       else
           return _render[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _render[k].set(v)
   end
}
setmetatable(render, render_mt)

local movie_mt = {
    __call = function (t,...)
        return makemovie(...)
    end,
   -- get method
   __index = function (t,k)
       if type(_render.movie[k]) == "function" then
           return _render.movie[k]
       else
           return _render.movie[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _render.movie[k].set(v)
   end
}
setmetatable(render.movie, movie_mt)

function renderstart(startframe, skipframe)
    render_startframe = startframe
    render_skipframe = skipframe
end

-- TODO: we may need to shift this down to the c_api in order to get the best
-- performance.
function rendermany(start, final, interval, ...)
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


function renderall(...)
    rendermany(0,smv.getfinalframe(),1,...)
end
