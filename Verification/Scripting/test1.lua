print("Running script for " .. fdsprefix .. ".")
--hidewindow()
print("Date: " .. os.date("%c"))
smv = require "smv"
ssf = require "ssf"
ssfparser = require "ssfparser"
string = require "string"

initsmvdata()
redWrite = function(...)
    io.stderr:write("\27[31m")
    io.stderr:write(...)
    io.stderr:write("\27[0m")
end
greenWrite = function(...)
    io.stderr:write("\27[32m")
    io.stderr:write(...)
    io.stderr:write("\27[0m")
end

tests = {run=0,passed=0,failed=0}

level = 0
function writeLevelIndent()
    for i=0,level,1 do
        io.stderr:write("  ")
    end
end

function display_test_results(tests)
    io.stderr:write("\tRun\tPassed\tFailed\n")
    io.stderr:write(string.format("Tests\t%d\t\27[32m%d\27[0m\t\27[31m%d\27[0m\n",
        tests.run, tests.passed, tests.failed))
end
-- test is function that must not throw an exception
-- if we won't to do some other test we simply place an assert within the
-- function being tested
function test (testName, test)
    writeLevelIndent()
    io.stderr:write(testName .. "\n")
    level = level + 1
    local flag,result = pcall(test)
    tests.run = tests.run + 1
    writeLevelIndent()
    if flag then
        tests.passed = tests.passed + 1
        io.stderr:write(":[")
        greenWrite("OK")
        io.stderr:write("]\n")
    else
        tests.failed = tests.failed + 1
        io.stderr:write(":[")
        redWrite("Failed")
        io.stderr:write("]\n")
        writeLevelIndent()
        io.stderr:write("Failed with an exception.\n")
        writeLevelIndent()
        io.stderr:write(result)
    end
    level = level - 1
end

-- test is function that must return an exception
function testException (testName, test)
    writeLevelIndent()
    io.stderr:write(testName .. "\n")
    level = level + 1
    local flag, result = pcall(test)
    tests.run = tests.run + 1
    writeLevelIndent()
    if not flag then
        tests.passed = tests.passed + 1
        io.stderr:write(":[")
        greenWrite("OK")
        io.stderr:write("]\n")
    else
        tests.failed = tests.failed + 1
        io.stderr:write(":[")
        redWrite("Failed")
        io.stderr:write("]\n")
        writeLevelIndent()
        io.stderr:write("This test should have throw an exception.\n")
    end
    level = level - 1
end
test("determine if script is running from within smokeview", function()
    -- determine whether the script is running from within Smokeview
    if smokeviewEmbedded then
        print("Embedded in smokeview")
    else
        print("Not embedded in smokeview")
    end
    assert(smokeviewEmbedded)
end)

testException("call error()", function()error()end)
test("print smokeview info", function()
    -- print information on smokeview
    io.write(string.format("Version: %s\n", smokeviewProgram.version))
    io.write(string.format("Git Hash: %s\n", smokeviewProgram.githash))
    io.write(string.format("Title: %s\n", smokeviewProgram.titlerelease))
    io.write(string.format("Build Date: %s\n", smokeviewProgram.builddate))
    io.write(string.format("FDS Git Hash: %s\n", smokeviewProgram.fdsgithash))
    io.write(string.format("Smokeview Path: %s\n", smokeviewProgram.smokeviewpath))
    io.write(string.format("Smokezip Path: %s\n", smokeviewProgram.smokezippath))
    io.write(string.format("Texture Directory: %s\n", smokeviewProgram.texturedir))
    
end)
-- print information on the model

-- load each different type of data file
test("load slice file", function() load.datafile("room_fire_02.sf") end)
-- disabled this test for speed
-- test("load boundary file", function() load.datafile("room_fire_01.bf") end)
-- test("load smoke3d file", function() load.datafile("room_fire_01.s3d") end)
testException("load compressed smoke3d file", function() load.datafile("room_fire_01.s3d.zs") end)
test("load particle file", function() load.datafile("room_fire.prt5") end)
testException("load non-existant file", function() load.datafile("abcdefg.hi") end)
test("load slice vector file", function() load.vdatafile("room_fire_01.sf") end)
testException("load boundary vector file", function() load.vdatafile("room_fire_01.bf") end)
testException("load smoke3d vector file", function() load.vdatafile("room_fire_01.s3d") end)
testException("load compressed smoke3d vector file", function() load.vdatafile("room_fire_01.s3d.sz") end)
testException("load particle vector file", function() load.vdatafile("room_fire_01.prt5") end)
testException("load non-existant vector file", function() load.vdatafile("qwert.yu") end)

-- unload all the loaded data
test("unload all data", function() unload.all()end)
test("window.size()", function() window.size(1024,768)end)
test("set render.dir", function() render.dir = "renders" end)

-- the following tests depend on data not being loaded
test("no loaded file tests", function()
    test("pre-reqs", function()unload.all()end)
    test("view.colorbar.flip1", function()
        test("set", function() view.colorbar.flip = true end)
        test("get", function() return view.colorbar.flip end)
        test("equal", function()
            assert(view.colorbar.flip == true, "get does not match set")
        end)
    end)
    -- test with another value in case default settings make in impact
    test("view.colorbar.flip2", function()
        test("set", function() view.colorbar.flip = false end)
        test("get", function() return view.colorbar.flip end)
        test("equal", function()
            assert(view.colorbar.flip == false, "get does not match set")
        end)
    end)
    test("view.colorbar.index1", function()
        test("set", function() view.colorbar.index = 128 end)
        test("get", function() return view.colorbar.index end)
        test("equal", function()
            assert(view.colorbar.index == 128, "get does not match set")
        end)
    end)
    test("view.colorbar.index2", function()
        test("set", function() view.colorbar.index = 197 end)
        test("get", function() return view.colorbar.index end)
        test("equal", function()
            assert(view.colorbar.index == 197, "get does not match set")
        end)
    end)
    -- test with an existing viewpoint
    test("view.viewpoint", function()
        local x = "internal"
        test("set", function() view.viewpoint = x end)
        tempyieldscript()
        test("get", function() return view.viewpoint end)
        test("equal", function()
            assert(view.viewpoint == x, "get does not match set\n get is: "
                   .. tostring(view.viewpoint) .. "\n  set is: " .. tostring(x))
        end)
    end)
    -- test with a non-existing viewpoint
    testException("view.viewpoint failure", function()
        view.viewpoint = "qwetmkl"
    end)
    test("window.width", function()
        test("set", function() window.width = 1234 end)
        test("get", function() return window.width end)
        test("equal", function()
            assert(window.width == 1234, "get does not match set")
        end)
    end)
    test("view.framenumber", function()
        local x = 33
        test("set", function() view.framenumber = x end)
        test("get", function() return view.framenumber end)
        test("equal", function()
            assert(view.framenumber == x, "get does not match set\n get is: "
                   .. tostring(view.framenumber) .. "\n  set is: " .. tostring(x))
        end)
    end)
    test("window.height", function()
        test("set", function() window.height = 567 end)
        test("get", function() return window.height end)
        test("equal", function()
            assert(window.height == 567, "get does not match set")
        end)
    end)
    test("timebar.visibility1", function()
        test("set", function() timebar.visibility = true end)
        test("get", function() return timebar.visibility end)
        test("equal", function()
            assert(timebar.visibility == true, "get does not match set\n get is: "
                   .. tostring(timebar.visibility) .. "\n  set is: " .. tostring(true))
        end)
    end)
    test("timebar.visibility2", function()
        test("set", function() timebar.visibility = false end)
        test("get", function() return timebar.visibility end)
        test("equal", function()
            assert(timebar.visibility == false, "get does not match set\n get is: "
                   .. tostring(timebar.visibility) .. "\n  set is: " .. tostring(false))
        end)
    end)
    -- test("timebar.visibility.toggle", function()
    --     local orig = timebar.visibility
    --     timebar.visibility.toggle()
    --     assert(timebat.visibility ~= orig)
    -- end)
    -- remember that this is being used with no data loaded
    -- also the time set is not guaranteed, so get may not always match set
    test("time", function()
        test("set", function() time = 127 end)
        test("get", function() return time end)
        test("equal", function()
            assert(time == 127, "get does not match set")
        end)
    end)
end)

    
-- the following tests depend on data being loaded
test("loaded file test", function()
    test("pre-reqs", function() load.datafile("room_fire_01.sf") end)
    test("load.tour", function() load.tour("Circular") end)
    test("tour.keyframe", function()
        test("set", function() tour.keyframe = 5 end)
        test("get", function() return tour.keyframe end)
        test("equal", function()
            assert(tour.keyframe == 5, "get does not match set")
        end)
    end)
    -- test("tour.view", function()
    --     test("set", function() tour.keyframe = 5 end)
    --     test("get", function() return tour.keyframe end)
    --     test("equal", function()
    --         assert(tour.keyframe == 5, "get does not match set")
    --     end)
    -- end)
    test("unload.tour", function() unload.tour() end)
    test("render", function() render() end)
    test("render.type get/set", function()
        test("set", function() render.type = "JPG" end)
        test("get", function() return render.type end)
        test("equal", function()
            assert(render.type == "JPG", "get does not match set: " .. render.type)
        end)
    end)
    testException("render.type invalid", function()
        render.type = "qwerr"
    end)
    test("render.movie", function()
        render.movie("testname", "testbase", 30)
    end)
    test("render.movie.type get/set", function()
        test("set", function() render.movie.type = "MP4" end)
        test("get", function() return render.movie.type end)
        test("equal", function()
            assert(render.movie.type == "MP4", "get does not match set")
        end)
    end)
    testException("render.movie.type invalid", function()
        render.movie.type = "qwer"
    end)
    test("model.chid", function() assert(model.chid) end)
    test("model.slices", function() assert(model.slices) end)
    test("number of slices", function()
        assert(model.nslices == 8, "number of slices is incorrect, is : "
            .. tostring(model.nslices))
    end)
    test("model.meshes", function() assert(model.meshes) end)
    test("number of meshes", function()
        assert(model.nmeshes == 1, "number of meshes is incorrect, is : "
            .. tostring(model.nmeshes))
    end)
    test("print slice info", function()
        print("printing slice info")
        print(model.slices)
        for key,value in pairs(model.slices) do print(key,value.label) end
        for key,value in pairs(model.slices) do print(key,value.file) end
    end)
end)
display_test_results(tests)
exit()
camera_mod_elev(45)
camera_mod_az(45)
camera_mod_az(-90)
settime(15.0)
-- yieldscript()
displayCB()
print("render14")
render()
camera_set_elev(90)
camera_set_az(0)
camera_set_projection_type(1)
camera_set_eyey(-3.078639)
render()
print("unloadall")
-- unloadall();
--loaddatavfile("room_fire_02.sf")
print("settime")
-- loaddatafile("room_fire_01.bf")
settime(20.0)
setframe(68)
render()
setframe(71)
settimebarvisibility(0)
camera_mod_eyex(0.1)
-- print("setrenderdir")
-- print("Current Script Directory: ")
--print(current_script_dir)
--print("Current Render Directory: ")
--print(current_render_dir)
-- yieldscript()
displayCB()
render()
print_times()
print("oh, here")
nframes = get_nglobal_times()
for i=0,nframes,1 do
    setframe(i)
    render()
end
print("nframes: " .. nframes)

initsliceinfo()
print(sliceinfo)
for key,value in pairs(sliceinfo) do print(key,value.label) end
print("Script for " .. fdsprefix .. " complete.")
exit()
