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
test("unload all data", function()unload.all()end)
test("window.size()", function()window.size(1024,768)end)
test("set render.dir", function()render.dir = "renders" end)

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
end)

display_test_results(tests)
    
exit()
-- the following tests depend on data being loaded
-- test("loaded file test", function()load.datafile("room_fire_01.sf"))
--hidewindow()
--loadinifile("room_fireX.ini")
setwindowsize(1024,768)
setrenderdir("renders")
--loadslice("TEMPERATURE", 1, 4.5)
loaddatafile("room_fire_01.sf")
--loaddatafile("room_fire_01.bf")
setcolorbarflip(0)
setcolorbarindex(128)
setgridvisibility(1)
settime(10.0)
-- yieldscript()
displayCB()
render()
--load3dsmoke("soot mass fraction")
setcolorbarflip(1)
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
