How to use the PyroGraph data processing utility.

** QUANTITIES OVERVIEW **
File name: quantities.csv
This csv formatted file stores information on the quantities that are used to generate scatter plot data. Each quantity represents a single scatter plot.  The plot legend labels (groups) and symbols are defined in the groups.csv and styles.csv files respectively.

* QUANTITIES PARAMETERS *
- Index, this is a unique numeric value Index that is used to refer to the quantity entry.  This quantity number will be referred to in the data sets configuration file.

- Quantity_Label, this is a shortened continuous character string used to identify the quantity.

- Scatter_Plot_Title, this is a LaTeX formatted string that is used for the plot title text.

- Ind_Title, this is the independent axis title in LaTeX format. 

- Dep_Title, this is the independent axis title in LaTeX format.

- Plot_Min, the minimum axis value for both x and y axis.

- Plot_Max, the maximum axis value for both x and y axis.

- Sigma_2_E, the calculated experimental error for bounding the quantity scatter data.

- Title_Position, where the title should be located on the plot. This is a pair of values that represent the percent of the x and y axis range from 0 to 1.  Ex. If the axis range is 0 to 200, and you want the title text to start at 10 x and 185 y the pair of values in this field would be [0.05,0.925] which is [(10/200),(185/200)].

- Key Position, acceptable locations are below; 
'upper right', 'upper left', 'lower left', 'lower right', 'right', 'center left', 'center right', 'lower center', 'upper center', 'center'

- Plot_Width, this is the width of the plot figure in inches. Only the width needs to be specified as scatter plots have a width to height ratio of 1:1, while comparison plots have a computed default 3:2 ratio.

- Plot_Filename, the path and filename for the resulting PDF scatter plot.  This path is appended to the 'output_directory' path variable that is set in the main script module.


** GROUPS OVERVIEW **
File name: groups.csv
This csv formatted file stores information on the groups that are used to organize scatter plot data.

* GROUPS PARAMETERS *
- Index, this is a unique numeric value Index that is used to refer to the group entry.  This group number will be referred to in the data sets configuration file.

- Group_Title, this is the character string that is used as the group label.  All points in a group on the scatter plot will be identified with this label in the legend and symbol style on the plot.

- Symbol_Style_Index, this is the 9 character Index code that is used to identify a style definition in the styles.csv file.  This links the group label to a style definition.


** STYLES OVERVIEW **
File name: styles.csv
This csv formatted file stores information on styles for plotting.
The plotter.py module will look to the parameters stored here to determine how to draw the lines and/or points on a plot.
NOTE: Scatter plots will never have lines drawn, so if a Style_Index is used for a scatter plot Group that has Line parameters defined, they will be ignored.  This allows reuse of styles for both scatter and comparison plots.

* STYLES PARAMETERS *
- Index, this is a unique numeric Index that is used to refer to a style definition.  This style number will be referred to in the groups configuration file, to associate a group with a particular plotting style.  It will also be used in the data sets configuration file to assign styles to comparison plot lines.

- Symbol_Style, defines the shape of the symbol to use if required. If it is a line only plot and no symbol is required, then place 'na' in this field.
The symbol style field accepts the following values shown in quotes;
'.' (point)
',' (pixel)
'o' (circle)
'v' (triangle_down)
'^' (triangle_up)
'<' (triangle_left)
'>' (triangle_right)
'1' (tri_down)
'2' (tri_up)
'3' (tri_left)
'4' (tri_right)
's' (square)
'p' (pentagon)
'*' (star)
'h' (hexagon1)
'H' (hexagon2)
'+' (plus)
'x' (x)
'D' (diamond)
'd' (thin_diamond)
'|' (vline)
'_' (hline)

- Edge_Color, If a symbol is used, and you want the symbol to have a different color for the edge than the fill color, or only want an edge to be colored and no fill is required, set the color for the symbol edges here.

- Fill_Color, this sets the fill color for symbols. For symbols that are only outlines and no fill color, set this field to 'na'.

- Symbol_Size, this is a value in points, that scales the size of the symbols. If a Symbol is not used, set this to 'na'.

- Line_Color, this sets the color for a line.  For Symbol only plots, set this to 'na', for line only plots, set the color here, and set Edge_Color and Fill_Color to 'na'.

- Line_Style, this sets the style of the line that is drawn.  If a line is not required, set this field to 'na'.
The line style field accepts the following styles; solid, dashed, dash_dot, dotted

- Line_Thickness, this is a value in points, that sets the thickness or width of the plotted line.  If a line is not plotted, this value should be 'na'.

