-- Run: texlua preview/test_serialization.lua
local repair
mode = "html-only"
Make = {
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
assert(repair('<math><mi mathvariant="bold">i</math></mi>') ==
       '<math><mi mathvariant="bold">i</mi></math>')
for _, gap in ipairs({'&#x00A0;', '&#xA0;', ' ', ' '}) do
  assert(repair('<a' .. gap .. 'id="page.527"></a>') == '<a id="page.527"></a>')
end
local valid = '<math><mo>&#x003E;</mo><mi>x</mi></math><a href="target">link</a>'
assert(repair(valid) == valid)
assert(repair(repair(valid)) == valid)
print('TeX4ht serialization regression checks passed.')
