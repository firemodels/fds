--- @module smv
local smv = {}
require "bounds"
require "clipping"
require "load"
require "render"
require "view"

-- set the defaults for renderall
render_startframe = 0
render_skipframe = 1

smv.getfinalframe = get_nglobal_times
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

return smv
