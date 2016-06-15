
local inioptions = {}

function length(table)
    local i = 0
    for key, val in pairs(table) do i = i + 1 end
    return i
end

function mockExec(name)
    print("Executing: ", name)
end

options =
    { AMBIENTLIGHT = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting ambient light to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.ambientlight = {r = r, g = g, b = b} end
        end
    , BACKGROUNDCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting background color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.backgroundcolor = {r = r, g = g, b = b} end
        end
    , BLOCKCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting block color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.blockcolor = {r = r, g = g, b = b} end
        end
    , BLOCKSHININESS = function(opt)
        local v = opt.argLines[1][1]
        print("Setting block shininess to " .. v)
        return function () view.color.blockshininess = v end
        end
    , BLOCKSPECULAR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting block specular to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.blockspecular = {r = r, g = g, b = b} end
        end
    , BOUNDCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting bound color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.boundcolor = {r = r, g = g, b = b} end
        end
    , COLORBAR = function(opt)
        local ncolors = opt.argLines[1][1]
        local texture_flag = opt.argLines[1][2]
        local contour_value = opt.argLines[1][3]
        local r, g, b
        local colors = {}
        for i=2,ncolors+1,1 do
            r = opt.argLines[i][1]
            g = opt.argLines[i][2]
            b = opt.argLines[i][3]
            table.insert(colors,{r=r,g=g,b=b})
        end
        return function ()
            -- print("COLORBAR not implemented")
            set_colorbar_colors(colors)
            -- view.colorbar.colors = colors
            set_colorbar_textureflag(textureflag ~= 0)
            -- view.colorbar.texture_flag = (textureflag ~= 0)
            setcolorbarindex(contour_value)
            -- view.colorbar.index = contour_value
            end
        end
    , COLORBAR_SPLIT = function(opt)
        return function ()
            print("COLORBAR_SPLIT not implemented")
            end
        end
    , GEOMSHOW = function(opt)
        local show_faces_interior = opt.argLines[1][1]
        local show_faces_exterior = opt.argLines[1][2]
        local show_faces_solid = opt.argLines[1][3]
        local show_faces_outline = opt.argLines[1][4]
        local smooth_geom_normal = opt.argLines[1][5]

        local show_volumes_interior = opt.argLines[2][1]
        local show_volumes_exterior = opt.argLines[2][2]
        local show_volumes_solid = opt.argLines[2][3]
        local show_volumes_outline = opt.argLines[2][4]

        local geom_vert_exag = opt.argLines[3][1]
        local geom_max_angle = opt.argLines[3][2]
        return function ()
            set_showfaces_interior(show_faces_interior)
            set_showfaces_exterior(show_faces_exterior)
            set_showfaces_solid(show_faces_solid)
            set_showfaces_outline(show_faces_outline)
            set_smoothgeomnormal(smooth_geom_normal)

            set_showvolumes_interior(show_volumes_interior)
            set_showvolumes_exterior(show_volumes_exterior)
            set_showvolumes_solid(show_volumes_solid)
            set_showvolumes_outline(show_volumes_outline)

            set_geomvertexag(geom_vert_exag)
            set_geommaxangle(geom_max_angle)
            end
        end
    , NORTHANGLE = function(opt)
        return function ()
            print("NORTHANGLE not implemented")
            end
        end
    , SHOWHRRLABEL = function(opt)
        return function ()
            print("SHOWHRRLABEL not implemented")
            end
        end
    , TREEPARMS = function(opt)
        return function ()
            print("TREEPARMS not implemented")
            end
        end
    , SCRIPTFILE = function(opt)
        return function ()
            print("SCRIPTFILE not implemented")
            end
        end
    , C_BOUNDARY = function(opt)
        return function ()
            print("C_BOUNDARY not implemented")
            end
        end
    , V_BOUNDARY = function(opt)
        return function ()
            print("V_BOUNDARY not implemented")
            end
        end
    , COLOR2BAR = function(opt)
        local ncolors = opt.argLines[1][1]
        local r, g, b
        local colors = {}
        for i=2,ncolors+1,1 do
            r = opt.argLines[i][1]
            g = opt.argLines[i][2]
            b = opt.argLines[i][3]
            table.insert(colors,{r=r,g=g,b=b})
        end
        return function ()
            view.color2bar.colors = colors
            end
        end
    , COLORBAR_FLIP = function(opt)
        local setting = opt.argLines[1][1]
        return function ()
            if (setting == 1) then
                setcolorbarflip(true)
            elseif (setting == 0) then
                setcolorbarflip(false)
            else
                error("invalid colorbar setting")
            end
            end
        end
    , DIFFUSELIGHT = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        print("Setting bound color to " .. r  .. " " .. g .. " "  .. b)
        return function () view.color.diffuselight = {r = r, g = g, b = b} end
        end
    , DIRECTIONCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_directioncolor(r,g,b) end
        end
    , FLIP = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_flip(setting) end
        end
    , FOREGROUNDCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_foregroundcolor(r,g,b) end
        end
    , HEATOFFCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_heatoffcolor(r,g,b) end
        end
    , HEATONCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_heatoncolor(r,g,b) end
        end
    , ISOCOLORS = function(opt)
        local shininess = opt.argLines[1][1]
        local default_opaqueness = opt.argLines[1][2]
        local specular = {}
        specular.r = opt.argLines[2][1]
        specular.g = opt.argLines[2][2]
        specular.b = opt.argLines[2][3]
        local nlevels = opt.argLines[3][1]
        local i
        local colors = {}
        for i=1,nlevels,1 do
            colors[i] = {r=opt.argLines[3+i][1],
                         g=opt.argLines[3+i][2],
                         b=opt.argLines[3+i][3]}
        end
        return function ()
            set_isocolors(shininess, default_opaqueness,
                          specular, colors)
        end
        end
    , COLORTABLE = function(opt)
        local ncolors = opt.argLines[1][1]
        local i
        local colors = {}
        for i=1,ncolors,1 do
            colors[opt.argLines[i+1][6]]
                 = {r=opt.argLines[1+i][1],
                         g=opt.argLines[1+i][2],
                         b=opt.argLines[1+i][3]}
        end
        return function ()
            set_colortable(colors)
        end
        end
    , LIGHT0 = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_light0(setting) end
        end
    , LIGHT1 = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_light1(setting) end
        end
    , LIGHTMODELLOCALVIEWER = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_light1(setting) end
        end
    , LIGHTMODELSEPARATESPECULARCOLOR = function(opt)
        local setting = opt.argLines[1][1]
        return function () set_lightmodelseparatespecularcolor(setting) end
        end
    , LIGHTPOS0 = function(opt)
        local x = opt.argLines[1][1]
        local y = opt.argLines[1][2]
        local z = opt.argLines[1][3]
        local w = opt.argLines[1][3]
        return function () set_lightpos0(x,y,z,w) end
        end
    , LIGHTPOS1 = function(opt)
        local x = opt.argLines[1][1]
        local y = opt.argLines[1][2]
        local z = opt.argLines[1][3]
        local w = opt.argLines[1][3]
        return function () set_lightpos1(x,y,z,w) end
        end
    , SENSORCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sensorcolor(r,g,b) end
        end
    , SENSORNORMCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sensornormcolor(r,g,b) end
        end
    , SETBW = function(opt)
        local geo_setting = opt.argLines[1][1]
        local data_setting = opt.argLines[1][2]
        return function () set_bw(geo_setting, data_setting) end
        end
    , SPRINKOFFCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sprinkleroffcolor(r,g,b) end
        end
    , SPRINKONCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_sprinkleroncolor(r,g,b) end
        end
    , STATICPARTCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_staticpartcolor(r,g,b) end
        end
    , TIMEBARCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_timebarcolor(r,g,b) end
        end
    , VENTCOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_ventcolor(r,g,b) end
        end
    , GRIDLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_gridlinewidth(v) end
        end
    , ISOLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_isolinewidth(v) end
        end
    , ISOPOINTSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_isopointsize(v) end
        end
    , LINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_linewidth(v) end
        end
    , PARTPOINTSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_partpointsize(v) end
        end
    , PLOT3DLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_plot3dlinewidth(v) end
        end
    , PLOT3DPOINTSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_plot3dpointsize(v) end
        end
    , SENSORABSSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_sensorabssize(v) end
        end
    , SENSORRELSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_sensorrelsize(v) end
        end
    , SLICEOFFSET = function(opt)
        local v = opt.argLines[1][1]
        return function () set_sliceoffset(v) end
        end
    , SMOOTHLINES = function(opt)
        local v = opt.argLines[1][1]
        return function () set_smoothlines(v) end
        end
    , SPHERESEGS = function(opt)
        local v = opt.argLines[1][1]
        return function () set_spheresegs(v) end
        end
    , SPRINKLERABSSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_sprinklerabssize(v) end
        end
    , STREAKLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_streaklinewidth(v) end
        end
    , TICKLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_ticklinewidth(v) end
        end
    , USENEWDRAWFACE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_usenewdrawface(v) end
        end
    , VECCONTOURS = function(opt)
        local v = opt.argLines[1][1]
        return function () set_veccontours(v) end
        end
    , VECLENGTH = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_veclength(a, b, c) end
        end
    , VECTORLINEWIDTH = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_vectorlinewidth(a, b, c) end
        end
    , VECTORPOINTSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_vectorpointsize(v) end
        end
    , VENTLINEWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_ventlinewidth(v) end
        end
    , VENTOFFSET = function(opt)
        local v = opt.argLines[1][1]
        return function () set_ventoffset(v) end
        end
    , WINDOWOFFSET = function(opt)
        local v = opt.argLines[1][1]
        return function () set_windowoffset(v) end
        end
    , WINDOWWIDTH = function(opt)
        local v = opt.argLines[1][1]
        return function () set_windowwidth(v) end
        end
    , WINDOWHEIGHT = function(opt)
        local v = opt.argLines[1][1]
        return function () set_windowheight(v) end
        end
    , BOUNDZIPSTEP = function(opt)
        local v = opt.argLines[1][1]
        return function () set_boundzipstep(v) end
        end
    , FED = function(opt)
        local v = opt.argLines[1][1]
        return function () set_fed(v) end
        end
    , FEDCOLORBAR = function(opt)
        local v = opt.argLines[1][1]
        return function () set_fedcolorbar(v) end
        end
    , ISOZIPSTEP = function(opt)
        local v = opt.argLines[1][1]
        return function () set_isozipstep(v) end
        end
    , NOPART = function(opt)
        local v = opt.argLines[1][1]
        return function () set_nopart(v) end
        end
    , PARTPOINTSIZE = function(opt)
        local v = opt.argLines[1][1]
        return function () set_partpointsize(v) end
        end
    , SHOWFEDAREA = function(opt)
        local v = opt.argLines[1][1]
        return function () set_showfedarea(v) end
        end
    , SLICEAVERAGE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_sliceaverage(a, b, c) end
        end
    , SLICEDATAOUT = function(opt)
        local a = opt.argLines[1][1]
        return function () set_slicedataout(a) end
        end
    , SLICEZIPSTEP = function(opt)
        local a = opt.argLines[1][1]
        return function () set_slicezipstep(a) end
        end
    , SMOKE3DZIPSTEP = function(opt)
        local a = opt.argLines[1][1]
        return function () set_smoke3dzipstep(a) end
        end
    , USER_ROTATE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[1][5]
        return function () set_userrotate(a, b, c, d, e) end
        end
    , APERTURE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_aperture(a) end
        end
    , AXISSMOOTH = function(opt)
        local a = opt.argLines[1][1]
        return function () set_axissmooth(a) end
        end
    , BLOCKLOCATION = function(opt)
        local a = opt.argLines[1][1]
        return function () set_blocklocation(a) end
        end
    , BOUNDARYTWOSIDE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_boundarytwoside(a) end
        end
    , CLIP = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_clip(a, b) end
        end
    , CONTOURTYPE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_contourtype(a, b) end
        end
    , CULLFACES = function(opt)
        local b = opt.argLines[1][2]
        return function () set_cullfaces(a) end
        end
    , ENABLETEXTURELIGHTING = function(opt)
        local a = opt.argLines[1][1]
        return function () set_texturelighting(a) end
        end
    , EYEVIEW = function(opt)
        local a = opt.argLines[1][1]
        return function () set_eyeview(a) end
        end
    , EYEX = function(opt)
        local a = opt.argLines[1][1]
        return function () set_eyex(a) end
        end
    , EYEY = function(opt)
        local a = opt.argLines[1][1]
        return function () set_eyey(a) end
        end
    , EYEZ = function(opt)
        local a = opt.argLines[1][1]
        return function () set_eyez(a) end
        end
    , FONTSIZE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_fontsize(a) end
        end
    , FRAMERATEVALUE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_frameratevalue(a) end
        end
    , GEOMDIAGS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_geomdiags(a, b, c) end
        end
    , GVERSION = function(opt)
        local a = opt.argLines[1][1]
        return function () set_gversion(a) end
        end
    , ISOTRAN2 = function(opt)
        local a = opt.argLines[1][1]
        return function () set_isotran2(a) end
        end
    , MESHVIS = function(opt)
        local n = opt.argLines[1][1]
        local i
        local vals = {}
        for i=1,n,1 do
            vals[i] = opt.argLines[1+i][1]
        end
        return function ()
            set_meshvis(vals)
        end
        end
    , OFFSETSLICE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_offsetslice(a) end
        end
    , OUTLINEMODE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_outlinemode(a, b) end
        end
    , P3DSURFACETYPE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_p3dsurfacetype(a) end
        end
    , P3DSURFACESMOOTH = function(opt)
        local a = opt.argLines[1][1]
        return function () set_p3dsurfacesmooth(a) end
        end
    , PROJECTION = function(opt)
        local a = opt.argLines[1][1]
        return function () set_projection(a) end
        end
    , SCALEDFONT = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[2][1]
        local e = opt.argLines[2][2]
        local f = opt.argLines[2][3]
        return function () set_scaledfont(a) end
        end
    , SHOWALLTEXTURES = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showalltextures(a) end
        end
    , SHOWAXISLABELS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showaxislabels(a) end
        end
    , SHOWBLOCKLABEL = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showblocklabel(a) end
        end
    , SHOWBLOCKS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showblocks(a) end
        end
    , SHOWCADANDGRID = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showcadandgrid(a) end
        end
    , SHOWCADOPAQUE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showcadopaque(a) end
        end
    , SHOWCEILING = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showceiling(a) end
        end
    , SHOWCOLORBARS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showcolorbars(a) end
        end
    , SHOWCVENTS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_showcvents(a) end
        end
    , SHOWDUMMYVENTS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showdummyvents(a) end
        end
    , SHOWEVACSLICES = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_showevacslices(a,b,c) end
        end
    , SHOWFLOOR = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showfloor(a) end
        end
    , SHOWFRAME = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showframe(a) end
        end
    , SHOWFRAMELABEL = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showframelabel(a) end
        end
    , SHOWFRAMERATE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showframerate(a) end
        end
    , SHOWFRAMERATE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showframe(a) end
        end
    , SHOWGRID = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showgrid(a) end
        end
    , SHOWGRIDLOC = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showgridloc(a) end
        end
    , SHOWHMSTIMELABEL = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showhmstimelabel(a) end
        end
    , SHOWHRRCUTOFF = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showhrrcutoff(a) end
        end
    , SHOWISO = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showiso(a) end
        end
    , SHOWISONORMALS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showisonormals(a) end
        end
    , SHOWLABELS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showlabels(a) end
        end
    , SHOWMEMLOAD = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showmemload(a) end
        end
    , SHOWOPENVENTS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_showopenvents(a, b) end
        end
    , SHOWOTHERVENTS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_showothervents(a, b) end
        end
    , SHOWSENSORS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_showsensors(a, b) end
        end
    , SHOWSLICEINOBST = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showsliceinobst(a) end
        end
    , SHOWSMOKEPART = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showsmokepart(a) end
        end
    , SHOWSPRINKPART = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showsprinkpart(a) end
        end
    , SHOWSTREAK = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        return function () set_showstreak(a,b,c,d) end
        end
    , SHOWTERRAIN = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showterrain(a) end
        end
    , SHOWTETRAS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_showtetras(a, b) end
        end
    , SHOWTHRESHOLD = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_showthreshold(a,b,c) end
        end
    , SHOWTICKS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showticks(a) end
        end
    , SHOWTIMEBAR = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtimebar(a) end
        end
    , SHOWTIMELABEL = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtimelabel(a) end
        end
    , SHOWTITLE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtitle(a) end
        end
    , SHOWTRACERSALWAYS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtracersalways(a) end
        end
    , SHOWTRIANGLES = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[1][5]
        local f = opt.argLines[1][6]
        return function () set_showtriangles(a,b,c,d,e,f) end
        end
    , SHOWTRANSPARENT = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtransparentvents(a) end
        end
    , SHOWTRANSPARENTVENTS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtransparentvents(a) end
        end
    , SHOWTRIANGLECOUNT = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtrianglecount(a) end
        end
    , SHOWVENTFLOW = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_showventflow(a,b,c) end
        end
    , SHOWVENTS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showvents(a) end
        end
    , SHOWWALLS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showwalls(a) end
        end
    , SKIPEMBEDSLICE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_skipembedslice(a) end
        end
    , SMOKESENSORS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_smokesensors(a,b) end
        end
    , STARTUPLANG = function(opt)
        local a = opt.argLines[1][1]
        return function () set_startuplang(a) end
        end
    , STEREO = function(opt)
        local a = opt.argLines[1][1]
        return function () set_stereo(a) end
        end
    , SURFINC = function(opt)
        local a = opt.argLines[1][1]
        return function () set_surfinc(a) end
        end
    , TERRAINPARMS = function(opt)
        local r_min = opt.argLines[1][1]
        local g_min = opt.argLines[1][2]
        local b_min = opt.argLines[1][3]
        local r_max = opt.argLines[2][1]
        local g_max = opt.argLines[2][2]
        local b_max = opt.argLines[2][3]
        local vert_factor = opt.argLines[3][1]
        return function () set_terrainparams(r_min, g_min, b_min, r_max, g_max,
                                             b_max, vert_factor) end
        end
    , TITLESAFE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_titlesafe(a) end
        end
    , TRAINERVIEW = function(opt)
        local a = opt.argLines[1][1]
        return function () set_trainerview(a) end
        end
    , TRANSPARENT = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_transparent(a, b) end
        end
    , TWOSIDEDVENTS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_twosidedvents(a, b) end
        end
    , VECTORSKIP = function(opt)
        local a = opt.argLines[1][1]
        return function () set_vectorskip(a) end
        end
    , VOLSMOKE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[1][5]
        local f = opt.argLines[2][1]
        local g = opt.argLines[2][1]
        local h = opt.argLines[2][1]
        local i = opt.argLines[2][1]
        local j = opt.argLines[2][1]
        local k = opt.argLines[2][1]
        local l = opt.argLines[2][1]
        return function () set_volsmoke(a, b, c, d, e, f, g, h, i, j, k, l) end
        end
    , ZOOM = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_zoom(a,b) end
        end
    , CELLCENTERTEXT = function(opt)
        local a = opt.argLines[1][1]
        return function () set_cellcentertext(a) end
        end
    , INPUT_FILE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_inputfile(a) end
        end
    , LABELSTARTUPVIEW = function(opt)
        local a = opt.argLines[1][1]
        return function () set_labelstartupview(a) end
        end
    , PIXELSKIP = function(opt)
        local a = opt.argLines[1][1]
        return function () set_pixelskip(a) end
        end
    , RENDERCLIP = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[1][5]
        return function () set_renderclip(a,b,c,d,e) end
        end
    , RENDERFILELABEL = function(opt)
        local a = opt.argLines[1][1]
        return function () set_renderfilelabel(a) end
        end
    , RENDERFILETYPE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_renderfiletype(a, b) end
        end
    , RENDEROPTION = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_renderfiletype(a, b) end
        end
    , UNITCLASSES = function(opt)
        local n = opt.argLines[1][1]
        local i
        local vals = {}
        for i=1,n,1 do
            -- table.insert(vals,i)
            vals[i] = opt.argLines[1+i][1]
        end
        return function ()
            set_unitclasses(vals)
        end
        end
    , ADJUSTALPHA = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_renderfiletype(a, b) end
        end
    , COLORBARTYPE = function(opt)
        local typ = opt.argLines[1][1]
        local label = opt.argLines[1][3]
        return function () set_colorbartype(typ, label) end
        end
    , EXTREMECOLORS = function(opt)
        local rmin = opt.argLines[1][1]
        local gmin = opt.argLines[1][2]
        local bmin = opt.argLines[1][3]
        local rmax = opt.argLines[1][4]
        local gmax = opt.argLines[1][5]
        local bmax = opt.argLines[1][6]
        return function ()
            set_extremecolors(rmin, gmin, bmin, rmax, gmax, bmax)
            end
        end
    , FIRECOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_firecolor(r,g,b) end
        end
    , FIRECOLORMAP = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_firecolormap(a,b) end
        end
    , FIREDEPTH = function(opt)
        local a = opt.argLines[1][1]
        return function () set_firedepth(a) end
        end
    , SHOWEXTREMEDATA = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_showextremedata(a,b,c) end
        end
    , SMOKECOLOR = function(opt)
        local r = opt.argLines[1][1]
        local g = opt.argLines[1][2]
        local b = opt.argLines[1][3]
        return function () set_smokecolor(r,g,b) end
        end
    , SMOKECULL = function(opt)
        local a = opt.argLines[1][1]
        return function () set_smokecull(a) end
        end
    , SMOKESKIP = function(opt)
        local a = opt.argLines[1][1]
        return function () set_smokeskip(a) end
        end
    , SMOKEALBEDO = function(opt)
        local a = opt.argLines[1][1]
        return function () set_smokealbedo(a) end
        end
    , SMOKERTHICK = function(opt)
        local a = opt.argLines[1][1]
        return function () set_smokerthick(a) end
        end
    , USEGPU = function(opt)
        local a = opt.argLines[1][1]
        return function () set_usegpu(a) end
        end
    , VOLSMOKE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[1][5]

        local f = opt.argLines[2][1]
        local g = opt.argLines[2][2]
        local h = opt.argLines[2][3]
        local i = opt.argLines[2][4]
        local j = opt.argLines[2][5]
        local k = opt.argLines[2][6]
        local l = opt.argLines[2][7]
        return function ()
            set_volsmoke(a,b,c,d,e,f,g,h,i,j,k,l)
            end
        end
    , SHOWHAZARDCOLORS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showhazardcolors(a) end
        end
    , SHOWHZONE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showhzone(a) end
        end
    , SHOWSZONE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showszone(a) end
        end
    , SHOWVZONE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showvzone(a) end
        end
    , SHOWZONEFIRE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showzonefire(a) end
        end
    , SHOWPATHNODES = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showpathnodes(a) end
        end
    , SHOWTOURROUTE = function(opt)
        local a = opt.argLines[1][1]
        return function () set_showtourroute(a) end
        end
    , TOURCOLORS = function(opt)
        local ar = opt.argLines[1][1]
        local ab = opt.argLines[1][2]
        local ac = opt.argLines[1][3]

        local br = opt.argLines[2][1]
        local bb = opt.argLines[2][2]
        local bc = opt.argLines[2][3]

        local cr = opt.argLines[3][1]
        local cb = opt.argLines[3][2]
        local cc = opt.argLines[3][3]

        local dr = opt.argLines[4][1]
        local db = opt.argLines[4][2]
        local dc = opt.argLines[4][3]

        local er = opt.argLines[5][1]
        local eb = opt.argLines[5][2]
        local ec = opt.argLines[5][3]

        local fr = opt.argLines[6][1]
        local fb = opt.argLines[6][2]
        local fc = opt.argLines[6][3]

        local gr = opt.argLines[7][1]
        local gb = opt.argLines[7][2]
        local gc = opt.argLines[7][3]
        return function ()
            set_tourcolors_selectedpathline(ar,ab,ac)
            set_tourcolors_selectedpathlineknots(br,bb,bc)
            set_tourcolors_selectedknot(cr,cb,cc)
            set_tourcolors_pathline(dr,db,dc)
            set_tourcolors_pathknots(er,eb,ec)
            set_tourcolors_text(fr,fb,fc)
            set_tourcolors_avatar(gr,gb,gc)
            end
        end
    , TOURCONSTANTVEL = function(opt)
        local a = opt.argLines[1][1]
        return function () set_tourconstantvel(a) end
        end
    , VIEWALLTOURS = function(opt)
        local a = opt.argLines[1][1]
        return function () set_viewalltours(a) end
        end
    , VIEWTIMES = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function () set_viewtimes(a,b,c) end
        end
    , VIEWTOURFROMPATH = function(opt)
        local a = opt.argLines[1][1]
        return function () set_viewtourfrompath(a) end
        end
    , AVATAREVAC = function(opt)
        local a = opt.argLines[1][1]
        return function () set_avatarevac(a) end
        end
    , GEOMETRYTEST = function(opt)
        local a = opt.argLines[1][1] -- int
        local b = opt.argLines[1][2] -- int
        local c = opt.argLines[1][3] -- float
        local d = opt.argLines[1][4] -- float

        local vals = {}
        for i=1,5,1 do
            vals[i] = opt.argLines[2][i]
        end
        for i=1,5,1 do
            vals[i+5] = opt.argLines[3][i]
        end

        local b1Vals = {}
        for i=1,6,1 do
            b1Vals[i] = opt.argLines[4][i]
        end

        local b2Vals = {}
        for i=1,6,1 do
            b2Vals[i] = opt.argLines[5][i]
        end
        for i=1,6,1 do
            b2Vals[i+6] = opt.argLines[6][i]
        end

        local b3Vals = {}
        for i=1,3,1 do
            b3Vals[i] = opt.argLines[7][i]
        end

        return function ()
            set_geometrytest(a, b, c, d, vals, b1Vals, b2Vals, b3Vals)
        end
        end
    , DEVICEVECTORDIMENSIONS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        return function () set_devicevectordimensions(a,b,c,d) end
        end
    , DEVICEBOUNDS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_devicebounds(a,b) end
        end
    , DEVICEORIENTATION = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        return function () set_deviceorientation(a,b) end
        end
    , GRIDPARMS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[2][1]
        local e = opt.argLines[2][2]
        local f = opt.argLines[2][3]
        return function ()
            set_gridparms(a,b,c,d,e,f)
            end
        end
    , GSLICEPARMS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[2][1]
        local f = opt.argLines[2][2]
        local g = opt.argLines[2][3]
        local h = opt.argLines[3][1]
        local i = opt.argLines[3][2]
        return function ()
            set_gridparms(a,b,c,d,{e,f,g},{h,i})
            end
        end
    , LOADFILESATSTARTUP = function(opt)
        local a = opt.argLines[1][1]
        return function ()
            set_loadfilesatstartup(a)
            end
        end
    , MSCALE = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        return function ()
            set_mscale(a,b,c)
            end
        end
    , SLICEAUTO = function(opt)
        local n = opt.argLines[1][1]
        local i
        local vals = {}
        for i=1,n,1 do
            vals[i] = opt.argLines[1+i][1]
        end
        return function ()
            set_sliceauto(vals)
        end
        end
    , MSLICEAUTO = function(opt)
        local n = opt.argLines[1][1]
        local i
        local vals = {}
        for i=1,n,1 do
            vals[i] = opt.argLines[1+i][1]
        end
        return function ()
            set_msliceauto(vals)
        end
        end
    , COMPRESSAUTO = function(opt)
        local a = opt.argLines[1][1]
        return function ()
            set_compressauto(a)
        end
        end
    , PART5PROPDISP = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[1][5]
        return function ()
            -- set_part5propdisp(a,b,c,d,e)
            print("PART5PROPDISP not implemented")
        end
        end
    , PART5COLOR = function(opt)
        return function ()
            print("PART5COLOR not implemented")
        end
        end
    -- , PART5CLASSVIS = function(opt)
    --     local n = opt.argLines[1][1]
    --     local i
    --     local vals = {}
    --     for i=1,n,1 do
    --         vals[i] = opt.argLines[1+i][1]
    --     end
    --     return function ()
    --         set_part5classvis(vals)
    --     end
    --     end
    , PROPINDEX = function(opt)
        local n = opt.argLines[1][1]
        local i
        local vals = {}
        for i=1,n,1 do
            vals[i] = {}
            vals[i][1] = opt.argLines[1+i][1]
            vals[i][2] = opt.argLines[1+i][2]
        end
        return function ()
            set_propindex(vals)
        end
        end
    , partclassdataVIS = function(opt)
        return function ()
            print("partclassdataVIS not implemented")
        end
        end
    , SHOOTER = function(opt)
        -- local n = opt.argLines[1][1]
        -- local i
        -- local vals = {}
        -- for i=1,n,1 do
        --     vals[i] = opt.argLines[1+i][1]
        -- end
        return function ()
            -- set_part5classvis(vals)
            print("SHOOTER not implemented")
        end
        end
-- SHOOTER
--  0.500000 0.000000 0.018868
--  0.250000 0.000000 0.000000
--  0.000000 0.000000 0.000000
--  1.000000 0.000000 4.000000
--  10 1 100 0 0
--  1.000000 1.000000
    , SHOWDEVICES = function(opt)
        local n = opt.argLines[1][1]
        local i
        local names = {}
        for i=1,n,1 do
            names[i] = opt.argLines[1+i][1]
        end
        return function ()
            set_showdevices(names)
        end
        end
    , SHOWDEVICEVALS = function(opt)
        local a = opt.argLines[1][1]
        local b = opt.argLines[1][2]
        local c = opt.argLines[1][3]
        local d = opt.argLines[1][4]
        local e = opt.argLines[1][5]
        local f = opt.argLines[1][6]
        local g = opt.argLines[1][7]
        local h = opt.argLines[1][8]
        return function ()
            set_showdevicevals(a,b,c,d,e,f,g,h)
        end
        end
    , SHOWMISSINGOBJECTS = function(opt)
        local a = opt.argLines[1][1]
        return function ()
            set_showmissingobjects(a)
        end
        end
    , TOURINDEX = function(opt)
        local a = opt.argLines[1][1]
        return function ()
            set_tourindex(vals)
        end
        end
    , USERTICKS = function(opt)
        -- local n = opt.argLines[1][1]
        -- local i
        -- local vals = {}
        -- for i=1,n,1 do
        --     vals[i] = opt.argLines[1+i][1]
        -- end
        return function ()
            -- set_part5classvis(vals)
            print("USERTICKS not implemented")
        end
        end
-- USERTICKS
--  0 1 5 1 1 1
--  0.000000 -1.000000 0.000000
--  0.000000 -1.000000 0.000000
--  424.000000 203.000000 16.000000
--  1.000000 1.000000 1.000000
--  1 1 1
    , XYZCLIP = function(opt)
        -- local n = opt.argLines[1][1]
        -- local i
        -- local vals = {}
        -- for i=1,n,1 do
        --     vals[i] = opt.argLines[1+i][1]
        -- end
        return function ()
            -- set_part5classvis(vals)
            print("XYZCLIP not implemented")
        end
        end
-- XYZCLIP
--  2
--  0 -0.424000 0 424.424011
--  0 94.497025 0 203.203995
--  0 -0.016000 1 5.193398

    , C_PARTICLES = function(opt)
        local minFlag  = opt.argLines[1][1]
        local minValue = opt.argLines[1][2]
        local maxFlag  = opt.argLines[1][3]
        local maxValue = opt.argLines[1][4]
        local label
        if (#opt.argLines[1] >= 5)
            then label = opt.argLines[1][5]
            else label = ""
            end
        return function ()
            set_c_particles(minFlag, minValue, maxFlag, maxValue, label)
        end
        end
    , C_PLOT3D = function(opt)
        -- local n = opt.argLines[1][1]
        -- local i
        -- local vals = {}
        -- for i=1,n,1 do
        --     vals[i] = opt.argLines[1+i][1]
        -- end
        return function ()
            -- set_part5classvis(vals)
            print("C_PLOT3D not implemented")
        end
        end
-- C_PLOT3D
--  5
--  1 0 1.000000 0 -0.000000
--  2 0 1.000000 0 -0.000000
--  3 0 1.000000 0 -0.000000
--  4 0 1.000000 0 -0.000000
--  5 0 1.000000 0 -0.000000
    , C_SLICE = function(opt)
        local minFlag  = opt.argLines[1][1]
        local minValue = opt.argLines[1][2]
        local maxFlag  = opt.argLines[1][3]
        local maxValue = opt.argLines[1][4]
        local label
        if (#opt.argLines[1] >= 5)
            then label = opt.argLines[1][5]
            else label = ""
            end
        return function ()
            set_c_slice(minFlag, minValue, maxFlag, maxValue, label)
        end
        end
    , CACHE_BOUNDARYDATA = function(opt)
        local v  = opt.argLines[1][1]
        return function ()
            set_cache_boundarydata(v)
        end
        end
    , CACHE_QDATA = function(opt)
        local v  = opt.argLines[1][1]
        return function ()
            set_cache_qdata(v)
        end
        end
    , PATCHDATAOUT = function(opt)
        local outputFlag  = opt.argLines[1][1]
        local tmin  = opt.argLines[1][2]
        local tmax  = opt.argLines[1][3]
        local xmin  = opt.argLines[1][4]
        local xmax  = opt.argLines[1][5]
        local ymin  = opt.argLines[1][6]
        local ymax  = opt.argLines[1][7]
        local zmin  = opt.argLines[1][8]
        local zmax  = opt.argLines[1][9]
        return function ()
            set_patchdataout(outputFlag, tmin, tmax, xmin, xmax, ymin, ymax,
                             zmin, zmax)
        end
        end
    , PERCENTILELEVEL = function(opt)
        local v  = opt.argLines[1][1]
        return function ()
            set_percentilelevel(v)
        end
        end
    , TIMEOFFSET = function(opt)
        local v  = opt.argLines[1][1]
        return function ()
            set_timeoffset(v)
        end
        end
    , TLOAD = function(opt)
        local beginFlag  = opt.argLines[1][1]
        local beginVal  = opt.argLines[1][2]
        local endFlag    = opt.argLines[1][3]
        local endVal   = opt.argLines[1][4]
        local skipFlag   = opt.argLines[1][5]
        local skipVal  = opt.argLines[1][6]
        return function ()
            set_tload(beginFlag, beginVal, endFlag, endValue, skipFlag,
                      skiValue)
        end
        end
    , V_PARTICLES = function(opt)
        return function ()
            print("V_PARTICLES not implemented")
        end
        end
-- V_PARTICLES
--  0 1.000000 0 0.000000
    , V5_PARTICLES = function(opt)
        return function ()
            print("V5_PARTICLES not implemented")
        end
        end
-- V5_PARTICLES
--  0 1.000000 0 0.000000 Uniform
    , V_PLOT3D = function(opt)
        return function ()
            print("V_PLOT3D not implemented")
        end
        end
-- V_PLOT3D
--  5
--  1 0 1.000000 0 1.000000
--  2 0 1.000000 0 1.000000
--  3 0 1.000000 0 1.000000
--  4 0 1.000000 0 1.000000
--  5 0 1.000000 0 1.000000
    , V_SLICE = function(opt)
        local minFlag  = opt.argLines[1][1]
        local minValue  = opt.argLines[1][2]
        local maxFlag  = opt.argLines[1][3]
        local maxValue  = opt.argLines[1][4]
        local label  = opt.argLines[1][5]
        -- local colon  = opt.argLines[1][6]
        local lineMin  = opt.argLines[1][7]
        local lineMax  = opt.argLines[1][8]
        local lineNum  = opt.argLines[1][9]
        return function ()
            set_v_slice(minFlag, minValue, maxFlag, maxValue,
                        label, lineMin, lineMax, lineNum)
        end
        end
    , V_TARGET = function(opt)
        return function ()
            print("V_TARGET not implemented")
        end
        end
-- V_TARGET
--  0 1.000000 0 0.000000
    , VIEWPOINT5 = function(opt)
        -- local v  = opt.argLines[1][1]
        return function ()
            print("VIEWPOINT5 not implemented")
            -- set_cache_qdata(v)
        end
        end
-- VIEWPOINT5
--  0 10 2
--  0.490669 -2.257067 0.018868 1.000000 -2
--  0.000000 0.000000 0.000000 1
--  0.500000 0.240566 0.018868
--  0.000000 90.000000
--  1.000000 0.000000 0.000000 0.000000
--  0.000000 1.000000 0.000000 0.000000
--  0.000000 0.000000 1.000000 0.000000
--  0.000000 0.000000 0.000000 1.000000
--  2 0 0 0 0 0 1
--  -0.424000 -1.204000 -0.016000 424.424011 203.203995 5.193398
--  topDown
    }

-- an ssfType is either a string or one of the three number types:
-- boolean, int, and float. A value of 0 1 one will be classed as a boolean by
-- this function, even though it may function as an int or a float and so on.
-- There is the patter boolean > int > float, so if a float is required it is
-- perfectly fine for a boolean to be given, as boolean is simple the most
-- restrictive type it satisfies.
function isInt(n) return (n % 1) == 0 end
function isBool(n) return n == 1 or n == 0 end
function ssfNumberType(n)
    if isBool(n) then return "boolean"
    elseif isInt(n) then return "integer"
    else return "float"
    end
end
function ssfType(arg)
    local t = type(arg)
    if t == "string" then return t
    elseif t == "number"
        then return ssfNumberType(arg)
        else error("invalid type: " .. t)
    end
end

function testArgType(specArgType, arg)
    local usedArgType = ssfType(arg)
    if specArgType == "string" then return usedArgType == "string"
    elseif specArgType == "float" then return usedArgType == "boolean"
        or usedArgType == "integer" or usedArgType == "float"
    elseif specArgType == "integer" then return usedArgType == "boolean"
        or usedArgType == "integer"
    else return usedArgType == "boolean"
    end
end

-- does not fault on too many arguments being given, they are simply ignored
function testArgTypes(matchedCommand, v)
    local argTypes = matchedCommand.argTypes
    local args = v.args
    for i, argType in ipairs(argTypes) do
        local arg = args[i]
        local optional
        if type(argType) == "table"
            then optional = argType.optional
                 argType = argType.type
            else optional = false
        end
        if arg
            then
                assert(testArgType(argType, arg), "Argument " .. i .. " is \""
                .. arg .. "\" but should be of the type " .. argType)
            else
                if not optional then
                    local expected = #argTypes
                    local given = #args
                    error(v.command .. ": Insufficient arguments."
                        .. " Expected " .. expected .. " but was given "
                        .. given)
                end
        end
    end
    return true
end

function inioptions.validate(v, execute)
    local matchedOption = options[v.name]
    if (matchedOption)
        then
            local execFunc = matchedOption(v)
            if execute then
                return execFunc()
            else
                return execFun
            end
        else
            print("Option " .. v.name .. " not recognised.")
            -- error("Option " .. v.name .. " not recognised.")
    end
end

return inioptions
