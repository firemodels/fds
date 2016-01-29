render = {}
function _render(...)
    local fstArg = select(1,...)
    if (...) == nil then renderC(...)
    elseif type(fstArg) == "string" then
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
    elseif type(fstArg) == "function" then renderC(fstArg())
     -- by only using the fstArg we strip any
                               -- extraneous arguments
    -- there is no argument and we stick with the default naming
    elseif type(fstArg) == "nil" then renderC(...)
    else error(type(fstArg) .. " is an invalid argument to render()")
    end
    -- TODO: return information about the file produced or something
end

render_mt = {
    __call = function (t,...)
        return _render(...)
    end
}

setmetatable(render, render_mt)

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
