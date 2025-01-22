# Matlab to Python conversion

## 1. Use an AI code converter

This step will get you 90% of the way there.  Then there usually a few extra steps to get a working script.

#### Use pandas.read_csv

I'd like to standardize on using Pandas for reading csv files.  This is not a hard rule, but there are many options and it would help if we are consistent.  You can give this instruction to your AI prior to your conversion.  You should end up with
```
import pandas as pd
```
at the top of your script and you read the csv file with, for example,
```
fds_data = pd.read_csv(fds_file, skiprows=2, header=None) # note: your header row may be different and remember Python is 0 based
```

## 2. Import `fdsplotlib`

You will use this library for the Matlab equivalent of `plot_style` and `addverstr`.  At the top of your script add the following:

```
import fdsplotlib
```

To make sure `fdsplotlib` is visible to all your python scripts add the following to your `~/.bashrc`
```
# add PYTHONPATH for FDS scripts
export PYTHONPATH=$FIREMODELS/fds/Utilities/Python:$PYTHONPATH
```

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
