--- @module smv
local smv = {}
require "bounds"
require "clipping"
require "load"
require "render"
require "view"
require "tour"

-- set the defaults for renderall
render_startframe = 0
render_skipframe = 1

-- predefined axes
X = 1
Y = 2
Z = 3

orthogonal = 1
perspective = 0

smv.getfinalframe = function()return get_nglobal_times()-1 end
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
model = {}
_model = {
    chid = {
        get = function()
            return fdsprefix
        end,
        set = function()
            error("model.chid is read-only")
        end
    },
    slices = {
        get = function()
            return sliceinfo
        end,
        set = function()
            error("model.slices is read-only")
        end,
        -- __len = function()
        --     error("len called slices")
        -- end
    },
    -- TODO: provide this by overriding the len operator
    nslices = {
        get = function()
            return #model.slices + 1
        end,
        set = function()
            error("model.nslices is read-only")
        end,
    },
    meshes = {
        get = function()
            -- this relies on initsmvdata being called first
            return meshinfo
        end,
        set = function()
            error("model.meshes is read-only")
        end
    },
    nmeshes = {
        get = function()
            return #model.meshes + 1
        end,
        set = function()
            error("model.nmeshes is read-only")
        end,
    }
}
local model_mt = {
   -- get method
   __index = function (t,k)
       if type(_model[k]) == "function" then
           return _model[k]
       else
           return _model[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _model[k].set(v)
   end
}
setmetatable(model, model_mt)

timebar = {}
_timebar = {
    visibility = {
        get = function()
            return get_timebar_visibility()
        end,
        set = function(v)
            return set_timebar_visibility(v)
        end,
        -- toggle = function ()
        --     timebar.visibility = not timebar.visibility
        -- end
    },
}
local timebar_mt = {
   -- get method
   __index = function (t,k)
       if type(_timebar[k]) == "function" then
           return _timebar[k]
       else
           return _timebar[k].get()
       end
   end,
   -- set method
   __newindex = function (t,k,v)
       _timebar[k].set(v)
   end
}
setmetatable(timebar, timebar_mt)

camera = {}
function camera.get()
    local camera = {
        rotationType = camera_get_rotation_type(),
        rotationIndex = camera_get_rotation_index(),
        eyePos = {
            x = camera_get_eyex(),
            y = camera_get_eyey(),
            z = camera_get_eyez()
        },
        zoom = camera_get_zoom(),
        viewDir = camera_get_viewdir(),
        zAngle = {
            az = camera_get_az(),
            elev = camera_get_elev()
        }
    }
    return camera
end

function camera.print(camera)
    io.write(string.format("rotationType: %d\n", camera.rotationType))
    -- io.write(string.format("rotationIndex: %d\n", camera.rotationIndex))
    io.write(string.format("viewId: %s\n", camera.viewId))
    io.write(string.format("eyePos: (%f,%f,%f)\n",
        camera.eyePos.x, camera.eyePos.y, camera.eyePos.x))
    io.write(string.format("zoom: %f\n", camera.zoom))
    io.write(string.format("viewAngle: %f\n", camera.viewAngle))
    io.write(string.format("directionAngle: %f\n", camera.directionAngle))
    io.write(string.format("elevationAngle: %f\n", camera.elevationAngle))
    io.write(string.format("projectionType: %d\n", camera.projectionType))
    io.write(string.format("viewDir: (%f,%f,%f)\n",
        camera.viewDir.x, camera.viewDir.y, camera.viewDir.z))
    io.write(string.format("zAngle: (%f,%f)\n",
        camera.zAngle.az, camera.zAngle.elev))
    -- io.write(string.format("transformMatrix: %d\n", camera.transformMatrix))
    local clipping = camera.clipping
    io.write(string.format("clipping: "))
    if clipping == nil or clipping.mode == nil then
        io.write(string.format("mode: 0\n"))
    else
        io.write(string.format("mode: %d\n", clipping.mode))
    end
    io.write(string.format("  x: ", camera.clipping))
    if clipping == nil or clipping.x == nil or clipping.x.min == nil then
        io.write(string.format("-inf", camera.clipping))
    else
        io.write(string.format("%f", clipping.x.min))
    end
    io.write(" to ")
    if clipping == nil or clipping.x == nil or clipping.x.max == nil then
        io.write(string.format("+inf", camera.clipping))
    else
        io.write(string.format("%f", clipping.x.max))
    end
    io.write("\n")
    io.write(string.format("  y: ", camera.clipping))
    if clipping == nil or clipping.y == nil or clipping.y.min == nil then
        io.write(string.format("-inf", camera.clipping))
    else
        io.write(string.format("%f", clipping.y.min))
    end
    io.write(" to ")
    if clipping == nil or clipping.y == nil or clipping.y.max == nil then
        io.write(string.format("+inf", camera.clipping))
    else
        io.write(string.format("%f", clipping.y.max))
    end
    io.write("\n")
    io.write(string.format("  z: ", camera.clipping))
    if clipping == nil or clipping.z == nil or clipping.z.min == nil then
        io.write(string.format("-inf", camera.clipping))
    else
        io.write(string.format("%f", clipping.z.min))
    end
    io.write(" to ")
    if clipping == nil or clipping.z == nil or clipping.z.max == nil then
        io.write(string.format("+inf", camera.clipping))
    else
        io.write(string.format("%f", clipping.z.max))
    end
    io.write("\n")
        -- transformMatrix = {
        --      1.000000 0.000000 0.000000 0.000000
        --      0.000000 1.000000 0.000000 0.000000
        --      0.000000 0.000000 1.000000 0.000000
        --      0.000000 0.000000 0.000000 1.000000
        -- },
end

function camera.set(camera)
    if camera == nil then
        error("camera.set: camera does not exist")
    end
    camera_set_eyex(camera.eyePos.x)
    camera_set_eyey(camera.eyePos.y)
    camera_set_eyez(camera.eyePos.z)
    camera_set_zoom(camera.zoom)
    camera_set_viewdir(camera.viewDir.x, camera.viewDir.y, camera.viewDir.z)
    camera_set_projection_type(camera.projectionType)
    camera_set_elev(camera.zAngle.elev)
    camera_set_az(camera.zAngle.az)
end
time = {}
function time.set(time)
    -- TODO: determine if the time is available
    return settime(time)
end


return smv
