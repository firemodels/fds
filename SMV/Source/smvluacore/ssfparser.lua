--- @module ssfparser
local ssfparser = {}

local ssfcommands = require "ssfcommands"
lpeg = require "lpeg"
lpeg.locale(lpeg)

-- only execute if running in smokeiew, else do a dry-run
if smokeviewEmbedded then execute = true else execute = false end

function buildInst(inst)
    local command = inst[1]
    local args = inst[2]
    command = {command = inst[1], args = inst[2]}
    printInstruction(command)
    command = ssfcommands.validate(command, execute)
    return command
end

whitespace = lpeg.space^0
integer = lpeg.S("-+")^-1 * lpeg.R("09")^1 / tonumber
float = lpeg.S("-+")^-1 * lpeg.R("09")^1 * lpeg.P(".")^0 * lpeg.R("09")^0
    / tonumber
stringChar = lpeg.alnum + lpeg.punct -- the chars that are valid in a string
string = lpeg.C(stringChar^1)
filepath = string
eol = lpeg.P("\n") + lpeg.P("\r\n")
command = lpeg.C(lpeg.alnum^1)
comment = lpeg.space^0 * lpeg.P("//")
    * lpeg.P(lpeg.alnum + lpeg.punct + lpeg.S(" \t"))^0
    * #eol
argument = float + integer + string + lpeg.S(" \t")^1 + comment +
    lpeg.P(eol * (#(lpeg.space) + eol * #comment))
arguments = lpeg.Ct(argument^0) * eol^-1
instruction = lpeg.Ct(command * arguments) / buildInst
script = lpeg.Ct(instruction^0)

function printInstruction(inst)
    print("INSTRUCTION")
    io.write("  Command:\t" .. inst.command .. "\n")
    io.write("  Arguments:\t[")
    for key,val in pairs(inst.args) do
        io.write(val .. ", ")
    end
    io.write("]\n")
end

function printScript(script)
    for key, val in pairs(script) do
        printInstruction(val)
    end
end

-- exampleInstruction = "SETTIMEVAL\n 0 223.0\n 5 \n astring"
-- exampleComment = "// test comment\n"
-- exampleInstructions = "SETTIMEVAL\n 0 223.0\n 5\n a string\nNEXTCOMMAND\n// test\nCOMMANDAFTER\n"
-- testCommand = command:match(exampleInstruction)
-- print("testCommand:", testCommand)
-- testInstruction = instruction:match(exampleInstruction)
-- print("testInstruction:", testInstruction)
-- if(testInstruction) then printInstruction(testInstruction) end
-- testScript = script:match(exampleInstructions)
-- print("testScript:", testScript)
-- if(testScript) then printScript(testScript) end
-- testComment = comment:match(exampleComment)
-- print("testComment:", testComment)
function parseSSF(filepath)
    local f = io.open(filepath, "r")
    local input = f:read("*all")
    local parsedScript = script:match(input)
	io.close()
end

runSSF = parseSSF

function test()
    parseSSF("test.ssf")
end

return ssfparser
