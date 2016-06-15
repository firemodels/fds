bounds = {}
bounds.slices = {}
-- bounds.slices["VIS_C0.9H0.1"] = {}
-- bounds.slices["VIS_C0.9H0.1"].x = {}
-- function bounds.slices["VIS_C0.9H0.1"].x.set(xMin, xMax)
function bounds.slices.set(name, min, max)
    if (min == nil)
        then
            set_slice_bound_min(name, false, 0)
        else
            set_slice_bound_min(name, true, min)
    end
    if (max == nil)
        then
            set_slice_bound_max(name, false, 0)
        else
            set_slice_bound_max(name, true, max)
    end
end

function bounds.slices.get(name)
    local min = get_slice_bound_min(name)
    local max = get_slice_bound_max(name)
    return {min = min, max = max}
end
