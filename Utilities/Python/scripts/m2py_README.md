# Matlab to Python conversion

## 1. Use an AI code converter

This step will get you 90% of the way there.  Then there usually a few extra steps to get a working script.

## 2. Import `fdsplotlib`

You will use this library for the Matlab equivalent of `plot_style` and `addverstr`.  At the top of your script add the following:

```
import fdsplotlib
```

For now, we will assume that you are running your scripts from the `fds/Utilities/Python/` directory where `fdsplotlib` lives.  If you need to run your python scripts from another directory, you need to add `fds/Utilities/Python/` to the path or else copy `fdsplotlibl.py` to your working directory.

## 3. Add the plot style parameters

Somewhere above where you make the plots add,

```
plot_style = fdsplotlib.get_plot_style("fds")
plt.rcParams["font.family"] = plot_style["Font_Name"]
plt.rcParams["font.size"] = plot_style["Label_Font_Size"]
# print(plot_style) # uncomment to see list of parameters
```

This will polulate a Python "dictionary" called `plot_style` with all the parameters.  To see a list of the parameters uncomment the print statement.

## 4. Add Git version string

Similar to the Matlab plots, we need to add

```
fdsplotlib.add_version_string(ax, git_file, plot_type='linear')
```

where `ax` is the handle of the plot.
