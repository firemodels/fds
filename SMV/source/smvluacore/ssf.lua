--- @module ssf
local ssf = {}

ssf.parser = require "ssfparser"

-- TODO: write parser for ssf files

--
-- function translate(v) =
--     if v.command  == "CBARFLIP" then return colorbarflip -- implemented
--     elseif v.command  == "CBARNORMAL") == 1)return SCRIPT_CBARNORMAL; -- implemented
--     elseif v.command  == "EXIT") == 1)return SCRIPT_EXIT; -- implemented
--     elseif v.command  == "KEYBOARD") == 1)return SCRIPT_KEYBOARD; -- TODO: requires a parser to interpret
--     elseif v.command  == "GSLICEORIEN") == 1)return SCRIPT_GSLICEORIEN; -- implemented
--     elseif v.command  == "GSLICEPOS") == 1)return SCRIPT_GSLICEPOS; -- implemented
--     elseif v.command  == "GSLICEVIEW") == 1)return SCRIPT_GSLICEVIEW; -- implemented
--     elseif v.command  == "LOAD3DSMOKE") == 1)return SCRIPT_LOAD3DSMOKE; -- implemented
--     elseif v.command  == "LOADBOUNDARY") == 1)return SCRIPT_LOADBOUNDARY; -- implemented
--     elseif v.command  == "LOADFILE") == 1)return SCRIPT_LOADFILE; -- implemented as loaddatafile
--     elseif v.command  == "LABEL") == 1)return SCRIPT_LABEL; -- implemented
--     elseif v.command  == "LOADINIFILE" return loaddatafile -- implemented
--     elseif v.command  == "LOADISO") == 1)return SCRIPT_LOADISO; -- implemented
--     elseif v.command  == "LOADPARTICLES") == 1)return SCRIPT_LOADPARTICLES; -- implemented
--     elseif v.command  == "LOADPLOT3D") == 1)return SCRIPT_LOADPLOT3D; -- implemented
--     elseif v.command  == "LOADSLICE") == 1)return SCRIPT_LOADSLICE; -- implemented
--     elseif v.command  == "LOADTOUR") == 1)return SCRIPT_LOADTOUR; -- implemented
--     elseif v.command  == "LOADVOLSMOKE") == 1)return SCRIPT_LOADVOLSMOKE; -- implemented
--     elseif v.command  == "LOADVOLSMOKEFRAME") == 1)return SCRIPT_LOADVOLSMOKEFRAME; -- implemented
--     elseif v.command  == "LOADVFILE") == 1)return SCRIPT_LOADVFILE; -- implemented
--     elseif v.command  == "LOADVSLICE") == 1)return SCRIPT_LOADVSLICE; -- implemented
--     elseif v.command  == "MAKEMOVIE") == 1)return SCRIPT_MAKEMOVIE; -- implemented (currently failing)
--     elseif v.command  == "PARTCLASSCOLOR") == 1)return SCRIPT_PARTCLASSCOLOR; -- implemented
--     elseif v.command  == "PARTCLASSTYPE") == 1)return SCRIPT_PARTCLASSTYPE; -- implemented
--     elseif v.command  == "PLOT3DPROPS") == 1)return SCRIPT_PLOT3DPROPS; -- implemented
--     elseif v.command  == "RENDERALL") == 1)return SCRIPT_RENDERALL; -- implemmented
--     elseif v.command  == "RENDERCLIP") == 1)return SCRIPT_RENDERCLIP; -- implemented
--     elseif v.command  == "RENDERDIR") == 1)return SCRIPT_RENDERDIR; -- implemented
--     elseif v.command  == "RENDERTYPE") == 1)return SCRIPT_RENDERTYPE; -- implemented
--     elseif v.command  == "MOVIETYPE") == 1)return SCRIPT_MOVIETYPE; -- implemented
--     elseif v.command  ==  "RENDERSIZE") == 1)return SCRIPT_RENDERSIZE; -- implemented
--     elseif v.command  == "RENDERDOUBLEONCE") == 1)return SCRIPT_RENDERDOUBLEONCE;
--     elseif v.command  == "RENDERONCE") == 1)return SCRIPT_RENDERONCE; -- implemented as render
--     elseif v.command  == "RENDERSTART") == 1)return SCRIPT_RENDERSTART; -- implemented
--     elseif v.command  == "SCENECLIP") == 1)return SCRIPT_SCENECLIP; -- implemented
--     elseif v.command  == "SETTOURKEYFRAME") == 1)return SCRIPT_SETTOURKEYFRAME; -- implemented
--     elseif v.command  == "SETTOURVIEW") == 1)return SCRIPT_SETTOURVIEW; -- implemented
--     elseif v.command  == "SETTIMEVAL") == 1)return SCRIPT_SETTIMEVAL; -- implemented
--     elseif v.command  == "SETVIEWPOINT") == 1)return SCRIPT_SETVIEWPOINT; -- implemented
--     elseif v.command  == "SHOWPLOT3DDATA") == 1)return SCRIPT_SHOWPLOT3DDATA; -- TODO: see TODO in C api
--     elseif v.command  == "UNLOADALL") == 1)return SCRIPT_UNLOADALL; -- implemented
--     elseif v.command  == "UNLOADTOUR") == 1)return SCRIPT_UNLOADTOUR; -- implemented
--     elseif v.command  == "VOLSMOKERENDERALL") == 1)return SCRIPT_VOLSMOKERENDERALL;
--     elseif v.command  ==  "ISORENDERALL")==1)return SCRIPT_ISORENDERALL;
--     elseif v.command  ==  "XSCENECLIP")==1)return SCRIPT_XSCENECLIP; -- implemented
--     elseif v.command  == "YSCENECLIP") == 1)return SCRIPT_YSCENECLIP;  -- implemented
--     elseif v.command  == "ZSCENECLIP") == 1)return SCRIPT_ZSCENECLIP; -- implemented
--     else return error
--     end
-- end
return ssf
