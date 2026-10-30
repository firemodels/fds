-- Run: texlua preview/test_serialization.lua
local repair
local last_interaction
mkutils = {file_in_builddir = function(path) return path end}
mode = "html-only"
Make = {
  testlogfile = function(parameters)
    last_interaction = parameters.interaction
    return 17
  end,
  add = function() end,
  tex4ht = function() end,
  match = function(self, pattern, callback)
    if not repair then repair = callback end
  end
}
package.preload["make4ht-filter"] = function()
  return function(filters) return filters[1] end
end
package.preload["make4ht-domfilter"] = function()
  return function() return function(content) return content end end
end
local directory = arg[0]:match("^(.*[/\\])") or ""
dofile(directory .. "full-guide.mk4")
local temporary = os.tmpname()
local parameters = {input = temporary, interaction = 'batchmode'}
for _, sample in ipairs({
  {text = string.rep('An ordinary log line without an error\n', 60000), expected = 'errorstopmode'},
  {text = '! Undefined control sequence.\nl.12 \\bad', expected = 'batchmode'},
  {text = './chapter.tex:12: Undefined control sequence.', expected = 'batchmode'}
}) do
  local file = assert(io.open(temporary .. '.log', 'w'))
  file:write(sample.text)
  file:close()
  assert(Make.testlogfile(parameters) == 17)
  assert(last_interaction == sample.expected)
  assert(parameters.interaction == 'batchmode')
end
os.remove(temporary .. '.log')
os.remove(temporary)
assert(repair('<math><mi mathvariant="bold">i</math></mi>') ==
       '<math><mi mathvariant="bold">i</mi></math>')
for _, gap in ipairs({'&#x00A0;', '&#xA0;', ' ', ' '}) do
  assert(repair('<a' .. gap .. 'id="page.527"></a>') == '<a id="page.527"></a>')
end
local valid = '<math><mo>&#x003E;</mo><mi>x</mi></math><a href="target">link</a>'
assert(repair(valid) == valid)
assert(repair(repair(valid)) == valid)
print('TeX4ht serialization regression checks passed.')
