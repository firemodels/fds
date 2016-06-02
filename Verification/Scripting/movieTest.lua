print("Running script for " .. fdsprefix .. ".")
--hidewindow()
print("Date: " .. os.date("%c"))
package.path=package.path .. ";" .. "../../SMV/Build/gnu_linux_64/?.lua"
smv = require "smv"
-- ssf = require "ssf"
-- ssfparser = require "ssfparser"
string = require "string"

-- this initsmvdata is necessary to bring some data into the Lua interpreter
-- from the model. This is included here rather than doing in the Smokeview
-- code to increase separation. This will likely be removed in future versions.
initsmvdata()

function mkMovie()
    -- this is a very quick hack of a sript as a demonstration
    -- depends on ffmpeg and ImageMagick
    -- it currently needs to save to disk in order to combine the two different
    -- data sets from Smokeview, however, it is possible to stream images
    -- directly to ffmpeg.
    local f = assert(io.popen("uname", r))
    local sys = assert(f:read('*a'))
    f:close()
    local sep
    print(sys)
    if string.find(sys, "Linux") then sep = "/"
    else sep = "\\"
    end
    os.execute("mkdir renders")
    os.execute("mkdir renders" .. sep .. "tempslice")
    os.execute("mkdir renders" .. sep .. "smoke")
    os.execute("mkdir renders" .. sep .. "combined")
    render.type = "PNG"
    unload.all()
    local cam = {
        rotationType = 1,
        eyePos = {x = 0.481481, y = -1.278639, z = 0.222222},
        zoom =  1.0,
        viewAngle = 0,
        directionAngle = 0,
        elevationAngle = 0,
        projectionType = 0,
        viewDir  = {x = 0.481481, y = 0.500000, z = 0.222222},
        zAngle = {az = -50, elev = 55},
        transformMatrix = nil,
        clipping = nil
    }
    camera.print(cam)
    camera.set(cam)
    timebar.visibility = true
    local movframes = 500
    -- step 1: load all the necessary data
    load.datafile("room_fire_01.sf")
    load.datafile("room_fire_01.s3d")
    -- step 2: hide all the data from view
    local movieHndl = assert(io.popen(string.format(
        "ffmpeg"
        .. " -y"
        .. " -r 15"
        .. " -i pipe:0"
        .. " renders/testMovie.mp4"), "wb"))
    for i=0,movframes,1 do
        setframe(i)
        -- step 3: show the temperature data
        show_slices_showall()
        -- step 4: render the temperature data
        io.stderr:write("rendering temperature data\n")
        render.dir = "renders/tempslice"
        render(function() return tostring(view.framenumber) end)
        -- step 5: hide the temperature data (via hide all)
        show_slices_hideall()
        -- step 6: show the smoke data
        show_smoke3d_showall()
        -- step 7: render the smoke data
        io.stderr:write("rendering smoke data\n")
        render.dir ="renders/smoke"
        render(function() return tostring(view.framenumber) end)
        -- step 8: hide the smoke data (via hide all)
        show_smoke3d_hideall()
        -- TODO: use multiple pipes to pipe simultaneous images
        io.stderr:write(string.format("combining frame %d\n",i))
        local imgHndl = assert(io.popen(string.format(
            "montage"
            .. " renders/tempslice/%d.png"
            .. " renders/smoke/%d.png"
            .. " -tile 2x1"
            .. " -geometry +0+0"
            .. " png:-"
            , i, i, i), "rb"))
        local comb = imgHndl:read('*a')
        local f = assert(io.open("testfile.png", "wb"))
        f:write(comb)
        movieHndl:write(comb)
        imgHndl:close()
    end
    movieHndl:close()
    -- os.execute(string.format(
    --     "ffmpeg"
    --     .. " -y"
    --     .. " -r 15"
    --     .. " -i renders/combined/%%d.png"
    --     .. " renders/testMovie.mp4"))
    -- -- rendermany(0,movframes,1,function() return tostring(view.framenumber) end)
    -- unload.all()
    -- print("done load")
    -- rendermany(0,movframes,1,function() return tostring(view.framenumber) end)
    -- io.stderr:write("rendering complete\n")
    -- io.stderr:write("combining frames\n")
end
mkMovie()
exit()
-- load.datafile("room_fire_01.sf")
-- load.datafile("room_fire_01.s3d")
-- show_smoke3d_hideall()
-- show_slices_hideall()
-- show_smoke3d_showall()