print("Running script for " .. fdsprefix .. ".")
--hidewindow()
print("Date: " .. os.date("%c"))
package.path=package.path .. ";" .. "../../SMV/Build/gnu_linux_64/?.lua"
smv = require "smv"
-- ssf = require "ssf"
iniparser = require "iniparser"
string = require "string"

-- this initsmvdata is necessary to bring some data into the Lua interpreter
-- from the model. This is included here rather than doing in the Smokeview
-- code to increase separation. This will likely be removed in future versions.
initsmvdata()

parseINI("test.ini")
exit()