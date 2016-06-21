--- @module ssf
local ssf = {}

ssf.parser = require "ssfparser"

-- TODO: the parser currently needs revision for the following items.
-- function translate(v) =
--     elseif v.command  == "KEYBOARD") == 1)return SCRIPT_KEYBOARD; -- TODO: requires a parser to interpret
--     elseif v.command  == "MAKEMOVIE") == 1)return SCRIPT_MAKEMOVIE; -- implemented (currently failing)
--     elseif v.command  == "RENDERDOUBLEONCE") == 1)return SCRIPT_RENDERDOUBLEONCE;
--     elseif v.command  == "RENDERONCE") == 1)return SCRIPT_RENDERONCE; -- implemented as render
--     elseif v.command  == "SHOWPLOT3DDATA") == 1)return SCRIPT_SHOWPLOT3DDATA; -- TODO: see TODO in C api
--     elseif v.command  == "VOLSMOKERENDERALL") == 1)return SCRIPT_VOLSMOKERENDERALL;
--     elseif v.command  ==  "ISORENDERALL")==1)return SCRIPT_ISORENDERALL;
--     else return error
--     end
-- end
return ssf
