smv = require "smv"
string = require "string"

blockage_view_method(3)
render("outline only raw")
view.blockages.method = "outline_only"
render("outline only interface")
blockage_outline_color(1)
render("blockageokage color outline raw")
blockage_outline_color(2)
render("foreground color outline raw")
exit()