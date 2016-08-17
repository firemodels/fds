print("Running script for " .. fdsprefix .. ".")
--hidewindow()
print("Date: " .. os.date("%c"))
package.path=package.path .. ";" .. "../../SMV/Build/gnu_linux_64/?.lua"
smv = require "smv"
-- ssf = require "ssf"
-- ssfparser = require "ssfparser"
string = require "string"

initsmvdata()

colors = {
    {r = 0.000000, g = 0.000000, b = 1.000000},
    {r = 0.000000, g = 0.359375, b = 1.000000},
    {r = 0.000000, g = 0.718750, b = 1.000000},
    {r = 0.000000, g = 1.000000, b = 0.776875},
    {r = 0.000000, g = 1.000000, b = 0.562500},
    {r = 0.000000, g = 1.000000, b = 0.203125},
    {r = 0.171875, g = 1.000000, b = 0.000000},
    {r = 0.531250, g = 1.000000, b = 0.000000},
    {r = 0.890625, g = 1.000000, b = 0.000000},
    {r = 1.000000, g = 0.746032, b = 0.000000},
    {r = 1.000000, g = 0.380952, b = 0.000000},
    {r = 1.000000, g = 0.000000, b = 0.000000}
}

set_colorbar_colors(12, colors)
set_colorbar_colors(12, colors)
view.colorbar.colors = colors
colorsE = get_colorbar_colors()
print(colorsE[4].b)
print(view.colorbar.colors[4].b)
-- view.colorbar.colors = colors
set_color2bar_colors(12, colors)
set_color2bar_colors(12, colors)
view.color2bar = colors
colorsE = get_color2bar_colors()
print(colorsE[4].b)
print(view.color2bar[4].b)
-- view.color2bar = colors
exit()