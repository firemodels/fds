tour = {}
_tour = {
    keyframe = {
        -- get = function ()
        --     return get_rendertype()
        -- end,
        set = function(v)
            return settourkeyframe(v)
        end
    },
    -- view = {
    --     -- get = function ()
    --     --     return get_rendertype()
    --     -- end,
    --     set = function(v)
    --         return settourview(v)
    --     end
    -- }
}
local tour_mt = {
    -- __call = function (t,...)
    --     return _renderF(...)
    -- end,
   -- get method
   __index = function (t,k)
       if type(_tour[k]) == "function" then
           return _tour[k]
       else
           return _tour[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _tour[k].set(v)
   end
}
setmetatable(tour, tour_mt)
