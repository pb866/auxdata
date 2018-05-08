__precompile__()

"""
# Module pyp

Load data and plot with PyPlot.


# Data structures

- PlotData


# Functions

- rd_data
- lineplot
- load_plotdata
- plot_data
- rename_DF
- calc_errors
- setup_axes
- set_cs
- plt_DataWithErrors
- redef_err
- set_log
"""
module pyp


##################
###  PREAMBLE  ###
##################

# Load packages
try using DataFrames
catch
  Pkg.add("DataFrame")
  using DataFrame
end
try using PyPlot
catch
  Pkg.add("PyPlot")
  using PyPlot
end
try using LaTeXStrings
catch
  Pkg.add("LaTeXStrings")
  using LaTeXStrings
end
try using Parameters
catch
  Pkg.add("Parameters")
  using Parameters
end

# Load further modules
if all(LOAD_PATH.!=Base.source_dir())  push!(LOAD_PATH,Base.source_dir())  end
using fhandle: test_file
using make_plots: sel_ls
# Export public functions
export rd_data,
       lineplot,
       load_PlotData,
       plot_data,
       PlotData

### NEW TYPES
@with_kw mutable struct PlotData
  x::Array{Number,1}
  y::Array{Number,1}
  xuerr::Union{Void,Array{Number,1}}=nothing
  xlerr::Union{Void,Array{Number,1}}=nothing
  yuerr::Union{Void,Array{Number,1}}=nothing
  ylerr::Union{Void,Array{Number,1}}=nothing
  label::Union{Void,AbstractString}=""
  marker::Union{AbstractString,Int64}="None"
  dashes::Union{Tuple{Int64,Int64},Array{Int64,1}}=Int64[]
  colour::Union{Void,AbstractString}=nothing
  lw::Number=1.4
end


##########################
###  PUBLIC FUNCTIONS  ###
##########################

"""
    rd_data(ifile; \*\*kwargs)


Read data from text file `ifile` in the following format:

    # Use '#' as first character for comments or data to be ignored.
    # Optional: specify column header names for DataFrame keywords
    # with keyword `jlheaders` as:
    # jlheaders: x1, y1, y2, ..., yn
    # (You may use whitespace, commas (,), semicolons (;), or pipes (|) as
    #  separators in the list above)
    <data separated by whitespace or any character/string>

If columns don't have the same length, it can be specified whether the first or last
columns will be filled with `NaN`s.

The function uses several keyword arguments (\*\*kwargs) for more freedom in the
file format or the selection of data.

### \*\*kwargs

- `ix` (`Int64`; default: `1`): Column index for column in `ifile` holding the x data
  (default column name in output DataFrame: `x`). If `ix` is set to `0`, no x column
  is assigned and only y columns are used in the DataFrame.
- `iy` (default: `0`): Column index/indices of y data columns in `ifile`. If `iy`
  is set to `0`, all columns starting at column 2 are assigned as y columns
  (default column name(s) in output DataFrame: `y1` ... `yn`).
  Columns can be specified using an integer for the selection of a single column,
  ranges (`<n(low)>:<n(high)>`), or arrays with special selections (`[n1, n2, ..., nn]`)
  where column order can be rearranged.
- `headers` (`Bool`; default: `false`): If headers is set to `true`, you need to specify
  column header names for all columns for the output DataFrame as described above
  using the keyword `jlheaders`. All columns will be read in and saved the same
  order as in `ifile`.
- `SF` (default value: `1` (no scaling)): You can optionally apply scaling factors
  to y data. If `SF` is an integer, the scaling factor will be applied to all
  y columns. You can apply scaling to each y column individually by providing an
  array of scaling factors of `length number of columns - number of x column`. If
  you only want to scale certain column(s), set the scaling factors for these columns
  and use `1` otherwise.
- `sep` (default: `whitespace`): You can specify any column separator with the
  keyword charactar `sep`. Separators can be any unicode character (even special
  characters such as `≠` or `α`) or string series of unicode characters
  (including whitespace).
- `fill` (`Int64`; default: `"last"`): If the column length of the input file varies,
  the `first` or `last` columns of the file are filled with `NaN`s according to
  the keyword. If you have a file with shorter columns to the right and the left,
  you either need to rearrange columns in the original data file or try to work
  with a specifically defined separator `sep`.
- `ncol` (`String`; default: `-1`): Defines the number of columns (x + y columns) in a file.
  If set to a negative number, the number of columns is derived from the `jlheaders`
  array or, if obsolete, from the first non-comment line of the file. You should
  only have to set the number of columns, if you have columns of different length
  with leading missing numbers.
- `skip_header` (`Int64`; default: `0`): Define how many lines to skip at the beginning
  of a file in case comments aren't used
- `skip_footer` (`Int64`; default: `0`): Define how many lines to skip at the end
  of a file in case comments aren't used
"""
function rd_data(ifile; ix::Int64=1, iy=0, headers::Bool=false, SF=1, sep::String="",
  fill::String="last", ncol::Int64=-1, skip_header::Int64=0, skip_footer::Int64=0)
  # Read input file
  ifile = test_file(ifile) # check existence of file
  lines = String[]; colnames = String[] # initialise arrays
  open(ifile,"r") do f
    # Read file
    lines = readlines(f)
    # Extract column header names, if specified
    if headers == true
      icol = find([contains(line,"jlheaders:") for line in lines])[1]
      colnames = split(lines[icol],"jlheaders:")[2]
      colnames = replace(colnames,r",|;|\|"," ")
      colnames = split(colnames)
      # Make sure, columns aren't rearranged
      if ix > 1  ix = 1  end
      iy = 0
    end
    # Skip first lines of a file, if skip_header is set to integer > 0
    deleteat!(lines,1:skip_header)
    # Skip last lines of a file, if skip_footer is set to integer > 0
    deleteat!(lines,1+length(lines)-skip_footer:length(lines))
    # Find and delete comment lines
    del=find(startswith.(lines,"#"))
    deleteat!(lines,del)
    # Find and delete empty lines
    del=find(lines.=="")
    deleteat!(lines,del)
  end

  # Set number of columns
  if ncol < 0 && length(colnames) > 0
    ncol = length(colnames)
  elseif ncol < 0
    ncol = length(split(lines[1]))
  end
  # Determine number of y columns for default case
  if iy == 0  && ix == 0
    iy = 1:ncol
  elseif iy == 0  && ix > 0
    iy = 2:ncol
  end

  # Initilise x and y data
  if ix > 0  x = Float64[]  end
  y = Matrix{Float64}(0, length(iy))
  # Loop over data lines
  for line in lines
    # Split into columns
    if sep == ""
      # Assume whitespace as default separator
      raw = split(line)
    else
      # Use separator, if specified
      raw = split(line,sep)
    end
    # Check number of current columns against maximum number of columns
    if length(raw) > ncol
      println("WARNING! Number of columns read in greater than defined number of columns.")
      println("The $fill $(length(raw)-ncol) columns are ignored.")
      if lowercase(fill[1]) == "l"
        raw = raw[1:ncol]
      else
        raw = raw[length(raw)-ncol+1:end]
      end
    end
    # Save current line to respective data arrays
    if ix > 0  push!(x,float(raw[ix]))  end
    ydat = transpose(float.(raw[1:end .!= ix]).*SF)
    ix == 0? nx = 0: nx = 1
    if length(ydat)<ncol && fill == "last"
      for i=length(ydat)+1:ncol-nx  ydat = hcat(ydat,NaN)  end
    elseif length(ydat)<ncol && fill == "first"
      for i=length(ydat)+1:ncol-nx  ydat = hcat(NaN,ydat)  end
    end
    y = vcat(y,ydat)
  end

  # Generate output DataFrame
  if ix == 0  output = DataFrame()
  else        output = DataFrame(x = x)
  end
  for i = 1:length(iy)
    output[Symbol("y$i")] = y[:,i]
  end
  # Rename headers, if names were specified
  if headers == true
    if length(output) != length(colnames)
      println("Warning! Length of column names not equal to number of columns.")
      println("Using standard names x for first column and y1...yn for remaining columns.")
    else
      for i = 1:length(colnames)
        rename!(output,names(output)[i],Symbol(colnames[i]))
      end
    end
  end

  # Return file data as DataFrame
  return output
end #function rd_data


"""
    lineplot(pdata, \*args; \*\*kwargs)
Plot `pdata` from a DataFrame to a pdf using PyPlot with the specifications
given by `*args` and `**kwargs`.


# \*args

- `xlab`: String for x axis label (default: `model time / hours`)
- `ylab`: String for y axis label (default: `concentration / mlc cm\$^{-3}\$ s\$^{-1}\$`)
- `ti`: String with plot title (default: empty `""`)

# \*\*kwargs

- `err`: Choose between options
  - `0` (default): no errors shown in plot
  - `1`: Symmetrical error – value in column 4 of input matrix indicates the
    ± range around the experimental value
  - `2`: Unsymmetrical error – value in columns 4 and 5 of the input matrix
    indicate the lower and upper error of the experimental value
  - `3`: Shade error within a percentage in the range
    `exp. value ± (exp. value × error)`
  - `4`: Shade error within a factor in the range
    `1/error × exp value … error × exp. value`
  - `5`: Shade error range with lower and upper bounds in columns 4 and 5, respectively
- `logscale`: set to `"x"`, `"y"` or `"xy"`, to set x axis, y axis or both axes to
  log scale, respectively
- `cs`: Set colour scheme from function sel_ls, 3 schemes available:
  `"line"` (default), `"source"`, and `"sink"`
- `lw`: Set line width of graphs
- `lt`: Set line style of graphs (different set of solid, dashed, dash-dotted,
  and dotted styles available)
- `lc`: Select single colors with integers, ranges or arrays with indexes of the
  colors from the colour scheme `cs`
- `mticks`: Use string `"on"` to show minor ticks and `"off"` suppress display
  of minor ticks
- `nmxt`: Set interval value of minor x ticks
- `nmyt`: Set interval value of minor y ticks
- `tsc`: Tuple with scaling factors for length of major (first tuple entry) and
  minor (second tuple entry) ticks in relation to their width
- `xlims`: Tuple with minimum and maximum x values; if you want to change only
  one value, use `nothing` for the value you want to keep
- `ylims`: Tuple with minimum and maximum y values; if you want to change only
  one value, use `nothing` for the value you want to keep
- `figsiz`: Tuple with (x, y) dimensions in inch for figure size
- `fntsiz`: Value of font size for tick labels and legend, axis labels are increased
  by 2
- `frw`: Set line width of plot border lines, grid lines, and tick width.
- `ti_offset`, `ax_offset`, `leg_offset`: Offsets for fontsizes of title, axes
  labels, and legend, respectively, from default sizes. Numbers
  (positive or negative) will be added to fontsizes.
- `legpos`: position of the legend; choose from the following options:
  - `"best"` or `0` (default)
  - `"upper right"` or `1`
  - `"upper left"` or `2`
  - `"lower left"` or `3`
  - `"lower right"` or `4`
  - `"right"` or `5`
  - `"center left"` or `6`
  - `"center right"` or `7`
  - `"lower center"` or `8`
  - `"upper center"` or `9`
  - `"center"` or `10`
- `legcol` (integer): number of legend columns (default: `1`)
- `SF`: scaling factor for y data, which can be included in the y legend, e.g. as
  `\\\$10^{-\$(log10(SF))} y\\\$`
"""
function lineplot(pdata, xlab::Union{String, LaTeXString}="model time / hours",
                  ylab::Union{String, LaTeXString}="concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
                  ti::Union{String, LaTeXString}="";  err::Int64=0, logscale::String="",
                  cs::String="line", lw::Number=1.4, lt=nothing, lc=nothing,
                  mticks::String="on", nmxt::Number=0, nmyt::Number=0,
                  xlims=nothing, ylims=nothing,
                  figsiz::Tuple{Number,Number}=(6,4), fntsiz::Number=12,
                  frw::Number=1, tsc::Tuple{Number,Number}=(4.5,2.5),
                  ti_offset::Number=4, ax_offset::Number=2, leg_offset::Number=0,
                  legpos="best", legcol::Int64=1, SF::Number=1)
  # Start plot
  fig, ax = subplots(figsize=figsiz)
  # Scale y data
  pdata[:,2] .*= SF
  if length(pdata[1,:]) ≥  4 && err ≠ 3 && err ≠ 4  pdata[:,4] .*= SF  end
  if length(pdata[1,:]) == 5 && err ≠ 3 && err ≠ 4  pdata[:,5] .*= SF  end
  # Add empty label column to input matrix, if missing
  if length(pdata[1,:]) == 2
    lab = String[]
    for i = 1:length(pdata[:,1])  push!(lab,"")  end
    pdata = hcat(pdata,lab)
  end

  # Set line colours and types
  if lc == nothing  lc = 1:length(pdata[:,1])  end
  if lt == nothing  lt = 1:length(pdata[:,1])  end
  lstyle = sel_ls(cs = cs, nc = lc, nt = lt)


  # Plot graphs
  for i = 1:length(pdata[:,1])
    ax[:plot](pdata[i,1], pdata[i,2], linewidth=lw, dashes=lstyle[2][i],
              color=lstyle[1][i], label=pdata[i,3])
    if err == 1
      ax[:fill_between](pdata[i,1], pdata[i,2].-pdata[i,4], pdata[i,2].+pdata[i,4],
        color=lstyle[1][i], alpha=0.2)
    elseif err == 2
      ax[:fill_between](pdata[i,1], pdata[i,2].-pdata[i,4], pdata[i,2].+pdata[i,5],
        color=lstyle[1][i], alpha=0.2)
    elseif err == 3
      ax[:fill_between](pdata[i,1], pdata[i,2].-(pdata[i,4].*pdata[i,2]),
        pdata[i,2].+(pdata[i,4].*pdata[i,2]), color=lstyle[1][i], alpha=0.2)
    elseif err == 4
      ax[:fill_between](pdata[i,1], (1./pdata[i,4]).*pdata[i,2],
        pdata[i,4].*pdata[i,2], color=lstyle[1][i], alpha=0.2)
    elseif err == 5
      ax[:fill_between](pdata[i,1], pdata[i,4], pdata[i,5],
        color=lstyle[1][i], alpha=0.2)
    end
  end

  # Set axes
  if logscale == "x"
    xmin = 10^floor(log10(minimum(minimum.(pdata[:,1]))))
    xmax = 10^ceil(log10(maximum(maximum.(pdata[:,1]))))
    xlim(xmin,xmax)
    ax[:set_xscale]("log")
  elseif logscale == "y"
    if err == 1
      ymin = 10^floor(log10(minimum(map(m->m[1],findmin.(
             max.(0,pdata[:,2].-pdata[:,4]))))))
      ymax = 10^ceil(log10(maximum(map(m->m[1],findmax.(
             pdata[:,2].+pdata[:,4])))))
    elseif err == 2
      ymin = 10^floor(log10(minimum(map(m->m[1],findmin.(
             max.(0,pdata[:,2].-pdata[:,4]))))))
      ymax = 10^ceil(log10(maximum(map(m->m[1],findmax.(
             pdata[:,2].+pdata[:,5])))))
    elseif err == 3
      ymin = 10^floor(log10(minimum(map(m->m[1],findmin.(
             max.(0,pdata[:,2].-(pdata[:,4].*pdata[:,2])))))))
      ymax = 10^ceil(log10(maximum(map(m->m[1],findmax.(
             pdata[:,2].+(pdata[:,4].*pdata[:,2]))))))
    elseif err == 4
      ymin = 10^floor(log10(minimum(map(m->m[1],findmin.(
             max.(0,(1./pdata[:,4]).*pdata[:,2]))))))
      ymax = 10^ceil(log10(maximum(map(m->m[1],findmax.(
             pdata[:,4].*pdata[:,2])))))
    elseif err == 5
      ymin = 10^floor(log10(minimum(map(m->m[1],findmin.(max.(0,pdata[:,4]))))))
      ymax = 10^ceil(log10(maximum(map(m->m[1],findmax.(pdata[:,5])))))
    else
      ymin = 10^floor(log10(minimum(map(m->m[1],findmin.(max.(0,pdata[:,2]))))))
      ymax = 10^ceil(log10(maximum(map(m->m[1],findmax.(pdata[:,2])))))
    end
    ylim(ymin, ymax)
    ax[:set_yscale]("log")
  elseif logscale == "xy"
    xmin = 10^floor(log10(minimum(minimum.(pdata[:,1]))))
    xmax = 10^ceil(log10(maximum(maximum.(pdata[:,1]))))
    xlim(xmin,xmax)
    ymin = 10^floor(log10(minimum(minimum.(pdata[:,2]))))
    ymax = 10^ceil(log10(maximum(maximum.(pdata[:,2]))))
    ylim(ymin, ymax)
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")
  #=
  else
    if xlims == nothing  xlims = (xlim()[1], xlim()[2])
    elseif xlims[1] == nothing  xlims = (xlim()[1], xlims[2])
    elseif xlims[2] == nothing  xlims = (xlims[1], xlim()[2])
    end
    if ylims == nothing  ylims = (ylim()[1], ylim()[2])
    elseif ylims[1] == nothing  ylims = (ylim()[1], ylims[2])
    elseif ylims[2] == nothing  ylims = (ylims[1], ylim()[2])
    end
    =#
  end
  ax[:set_xlim](xlims); ax[:set_ylim](ylims);

  # Set plot title
  ax[:set_title](ti, fontsize=fntsiz+ti_offset)
  # Generate axes labels and legend
  ax[:set_xlabel](xlab,fontsize=fntsiz+ax_offset)
  ax[:set_ylabel](ylab,fontsize=fntsiz+ax_offset)
  if any(pdata[:,3].!="")
    ax[:legend](fontsize=fntsiz+leg_offset, loc=legpos, ncol=legcol)
  end

  # Format plot
  if mticks == "on"
    minorticks_on()
  else
    minorticks_off()
  end
  # Set minor x ticks
  if nmxt > 0
    mx = matplotlib[:ticker][:MultipleLocator](nmxt)
    ax[:xaxis][:set_minor_locator](mx)
  end
  # Set minor y ticks
  if nmyt > 0
    my = matplotlib[:ticker][:MultipleLocator](nmyt)
    ax[:yaxis][:set_minor_locator](my)
  end
  # Format ticks and frame
  Mtlen = tsc[1]⋅frw
  mtlen = tsc[2]⋅frw
  ax[:tick_params]("both", which="both", direction="in", top="on", right="on",
    labelsize=fntsiz, width=frw)
  ax[:tick_params]("both", which="major", length=Mtlen)
  ax[:tick_params]("both", which="minor", length=mtlen)
  ax[:grid](linestyle=":", linewidth = frw)
  ax[:spines]["bottom"][:set_linewidth](frw)
  ax[:spines]["top"][:set_linewidth](frw)
  ax[:spines]["left"][:set_linewidth](frw)
  ax[:spines]["right"][:set_linewidth](frw)
  tight_layout()

  # Reset y data
  pdata[:,2] ./= SF
  if length(pdata[1,:]) ≥  4 && err ≠ 3 && err ≠ 4  pdata[:,4] ./= SF  end
  if length(pdata[1,:]) == 5 && err ≠ 3 && err ≠ 4  pdata[:,5] ./= SF  end

  return fig
end #function lineplot


"""
    load_PlotData(pltdata::DataFrames.DataFrame;  \*\*kwargs)

Pack x and y data with possible assoiciated errors from DataFrame `pltdata`
as well as formatting parameters into a new DataType `PlotData`.


### \*\*kwargs

+ `err` (`String`):
  - `"None"` (no errors, **default**)
  - `"rangex"`, `"rangey"`, `"range"` (± x/y/x and y)
  - `"pmrangex"`, `"pmrangey"`, `"pmrange"` (± x/y/x and y with different lower/upper values)
  - `"percentx"`, `"percenty"`, `"percent"` (± err⋅x/y/x and y)
  - `"pmpercentx"`, `"pmpercenty"`, `"pmpercent"` (± err⋅x/y/x and y with different lower/upper values)
  - `"factorx"`, `"factory"`, `"factor"` (x/y/x and y ⋅1/err and ⋅err, respectively)
  - `"pmfactorx"`, `"pmfactory"`, `"pmfactor"` (as above with different lower/upper values)
  - `"valuex"`, `"valuey"`, `"value"` (err value directly taken from column)
+ `mt` (`Union{AbstractString,Int64}`): PyPlot marker type, see e.g. https://matplotlib.org/api/markers_api.html
  - `"None"` (**default**) to suppress use of markers
+ `lt` (`Union{String,Tuple{Int64,Int64},Array{Int64,1}}`), tuple or array with
  on/off PyPlot dash definitions
  - empty array means a solid line (**default**)
  - Tuple with 2 entries or array with 4 entries where odd entries are `0` suppresses lines
+ `lc` (`Union{Void,AbstractString}`) PyPlot line/marker colours (**default:** `nothing` –
  use PyPlot default colours), see e.g. https://matplotlib.org/examples/color/named_colors.html
+ `lw` (`Number`): linewidth (**default:** `1.4`)
+ `SF` (`Number`): scaling factor of y data and associated errors (**default:** `1`, no scaling)
+ `label` (`String`): Label for legend
  - `""` (empty string, **default**): no legend label
+ `renameDF` (`Union{Bool,Array{Symbol,1}}`):
  - `"true"` (assume columns in the order: `x`, `y`, `ylerr`, `yuerr`, `xlerr`, `xuerr`)
  - `"false"` (columns already in correct order with names as for `true`)
  - `Array{Symbol, 1}`: give array with column names for `x`, `y`, `ylerr`, `yuerr`, `xlerr`, `xuerr`
   (give complete list, even if columns are incomplete due to the choice of `err`)
"""
function load_PlotData(plotdata::DataFrames.DataFrame;  err::String="None",
         mt::Union{AbstractString,Int64}="None",
         lt::Union{String,Tuple{Int64,Int64},Array{Int64,1}}=Int64[],
         lc::Union{Void,AbstractString}=nothing, lw::Number=1.4, SF::Number=1,
         label::String="", renameDF::Union{Bool,Array{Symbol,1}}=true)

  # Make copy of plotdata that can be altered
  pltdata = deepcopy(plotdata)
  # (Re-)define column names of DataFrame
  if renameDF == true
    DFnames = Symbol[:x, :y, :ylerr, :yuerr, :xlerr, :xuerr]
    pltdata = rename_DF(pltdata, DFnames, err)
  elseif isa(renameDF,Array{Symbol, 1})
    if length(renameDF) ≠ 6
      println("\'renameDF\' not correctly defined. Define all column names for")
      println("x, y, ylerr, yuerr, xlerr, and xuerr. Script stopped."); exit()
    end
    DFnames = deepcopy(renameDF)
  end

  # Calculate error columns depending on choice of `err`
  pltdata = calc_errors(pltdata, DFnames, err)

  # Scale data
  for s in DFnames[2:4]  try pltdata[s] .*= SF end  end

  # Save errors
  errors = []
  if err ≠ "None"
    for s in DFnames[3:end]
      try push!(errors,pltdata[s])
      catch
        push!(errors,nothing)
      end
    end
  else
    errors = [nothing, nothing, nothing, nothing]
  end

  # Return PlotData type
  return PlotData(x = pltdata[DFnames[1]], y = pltdata[DFnames[2]],
         xuerr = errors[4], xlerr = errors[3], yuerr = errors[2], ylerr = errors[1],
         label = label, marker = mt, dashes = lt, colour = lc, lw = lw)

end #function load_PlotData


"""
    plot_data(plot_list::Union{Array{PlotData, 1}, Array{pyp.PlotData, 1}}, \*args, \*\*kwargs)

Generate scatter and/or line plots from the array `plot_data` of Type `PlotData`.


### \*args

+ `xlab` (`Union{String, LaTeXString}`): x axis label
  - **default:** `"model time / hours"`
+ `ylab` (`Union{String, LaTeXString, Array{String,1}, Array{LaTeXString,1}}`):
  y axis label, `Array` can be used, if a second axis with a different label is introduced
  - **default:** `"concentration / mlc cm\$^{-3}\$ s\$^{-1}\$"`
+ `ti` (`Union{String, LaTeXString}`): Plot title
  - **default:** `""` (empty string) for no title


### \*\*kwargs

For all arguments concerning y data the following applies:
- The dataset can be split into 2 subsets and a 2. y-axis with a different scale
  can be assigned (see optional parameter `twinax`)
- If a second axis is defined, parameters concerning y-data can be compiled in an
  array of length `2` with different values for the first and second y-axis,
  respectively
- If the second y-axis is defined and only a single parameter concerning y-data
  is defined, then this parameter applies to both axes

Keyword arguments for `function plot_data` are:
+ `twinax` (`Array{Int64,1}`): optional array to devide data into 2 subsets and
  assign second dataset to a 2. y-axis; use array of `length` of the x and y data
  with integers `1` and `2` to assign each data to axis 1 and 2
  (**default:**: empty array `Int64[]`, no 2. y-axis)
+ `logscale` (`Union{String, Array{String, 1}}`): use `"x"`, `"y"` or `"xy"` to
  set either x-, y- or both axes to a logarithmic scale;
  if 2. y-axis is present, arrays may be used to have different scales on each y-axis;
  (**default:** `""` – linear scale)
+ `cs` (`Union{String, Array{String,1}}`): Define colour scheme `line`, `source`
  or `sink` from `function sel_ls`; array with 2 different colour schemes may be
  used in case of 2. y-axis
  - `""` (**default:**): colours as defined in the `PlotData`
  or default PyPlot colour scheme)
  - `"own"`: set your own colour schemes using parameters `lt` and `lc`
   (mandatory, if `cs` is set to `"own"`)
+ `lt`: dash type (defined by PyPlot's `dashes`)
  - for use only, when `cs` is switched to `"own"`
  - Array of `length 2` with Arrays of `length(plot_data)` where inner arrays
   specify line type of each `PlotData`
  - use tuple or array of integers to define on/off dash pixels
  - `[]` (empty array, **default**): solid lines
  - tuple with 2 entries or array with at least 4 entries where odd entries are `0`
    for no lines
+ `lc`: line/marker colour (defined by PyPlot's `color` keywords, see e.g.
  https://matplotlib.org/examples/color/named_colors.html)
  - for use only, when `cs` is switched to `"own"`
  - Array of `length 2` with Arrays of `length(plot_data)` where inner arrays
   specify the line/marker colour of each `PlotData`
+ `mticks` (`String`): Switch minor ticks on/off (**default:** "on")
+ `nmxt` (`Int64`): Number of minor x ticks between major ticks
  (**default:** `0` – automatic determination by PyPlot)
+ `nmyt` (`Union{Int64, Array{Int64,1}}`): Number of minor y ticks between major ticks
  - `0` (**default**): automatic determination by PyPlot
  - If 2. y axis is defined, an array can be used to used different specifications
   for each axis
+ `Mxtint`/`Mytint` (`Number`/`Union{Number,Array{Number,1}}`; **default:** `0`)
  interval size of major ticks in x- and y-axis, respectively; asign different intervals
  to 2nd y-axis using an array of integers
+ `xlims`/`ylims`: Tuple of minimum/maximum value for each axis
  - `nothing` may be used for default PyPlot values (**default**)
  - `nothing` may be used for one value within the tuple to use the default,
   e.g. `ylims = (0, nothing)`
  - if 2. y-axis is defined, `ylims` may be an array of tuples with different
   specifications
+ `figsiz` (`Tuple{Number,Number}`): figure size width × height in inches
  (**default:** `(6,4)`)
+ `fntsiz` (`Number`): default font size (**default:** `12`)
+ `frw` (`Number`): default line width of outer axes frame (**default:** `1`)
+ `tsc` (`Tuple{Number,Number}`): Tuple with scaling factors for length of major
  (first tuple entry) and minor (second tuple entry) ticks in relation to their
  width (**default:** `(4.5,2.5)`)
+ `cap_offset` (`Number`): deviation from default cap size of `3` of marker error bars
  (**default:** `0`)
+ `ti_offset`, `ax_offset`, `leg_offset`: Offsets for fontsizes of title, axes
  labels, and legend, respectively, from default sizes. Numbers
  (positive or negative) will be added to fontsizes
  (**defaults:** `4`, `2`, `0`)
+ `axcol` (`Union{String,Array{String,1}}`): colour of y axis label and tick numbers;
  if 2. y axis is defined, different values may be defined in an array for each axis
  (**default:** `"black"`)
+ `legpos` (`Union{String, Int64, Array{String,1}, Array{Int64,1}}`):
  position of the legend; choose from the following options:
  - `"best"` or `0` (default)
  - `"upper right"` or `1`
  - `"upper left"` or `2`
  - `"lower left"` or `3`
  - `"lower right"` or `4`
  - `"right"` or `5`
  - `"center left"` or `6`
  - `"center right"` or `7`
  - `"lower center"` or `8`
  - `"upper center"` or `9`
  - `"center"` or `10`
- `legcol` (`Union{Int64, Array{Int64,1}}`): number of legend columns
  (**default:** `1`)


### To Do:

- Optional marker type (add to sel_ls)
"""
function plot_data(plot_list::Union{Array{PlotData, 1}, Array{pyp.PlotData, 1}},
                  xlab::Union{String, LaTeXString}="model time / hours",
                  ylab::Union{String, LaTeXString, Array{String,1}, Array{LaTeXString,1}}=
                  "concentration / mlc cm\$^{-3}\$ s\$^{-1}\$",
                  ti::Union{String, LaTeXString}=""; twinax::Array{Int64,1}=Int64[],
                  logscale::Union{String, Array{String, 1}}="",
                  cs::Union{String, Array{String,1}}="", lt=[], lc="black",
                  mticks::String="on", nmxt::Int64=0, nmyt::Union{Int64, Array{Int64,1}}=0,
                  Mxtint::Number=0, Mytint::Union{Number,Array{Number,1}}=0,
                  xlims=nothing, ylims=nothing,
                  figsiz::Tuple{Number,Number}=(6,4), fntsiz::Number=12,
                  frw::Number=1, tsc::Tuple{Number,Number}=(4.5,2.5), cap_offset::Number=0,
                  ti_offset::Number=4, ax_offset::Number=2,
                  axcol::Union{String,Array{String,1}}="black", leg_offset::Number=0,
                  legpos::Union{String, Int64, Array{String,1}, Array{Int64,1}}="best",
                  legcol::Union{Int64, Array{Int64,1}}=1)

  # Start plot
  fig, ax1 = subplots(figsize=figsiz)
  # Check for twin axes and devide datasets
  plt, ax2, ax_2, ylab, logscale, xlims, ylims, Mytint, nmyt, cs, axcol, legpos, legcol =
    setup_axes(plot_list, twinax, ylab, logscale, xlims, ylims, Mytint, nmyt,
    cs, axcol, legpos, legcol)
  if ax2 ≠ nothing  ax = [ax1, ax2]
  else              ax = [ax1]
  end

  # set colour scheme
  plt = set_cs(plt, cs, lc, lt)

  # Plot data and associated errors
  plt[1], ax1 = plt_DataWithErrors(plt[1], ax1, cap_offset)
  if ax_2  plt[2], ax2 = plt_DataWithErrors(plt[2], ax2, cap_offset)  end

  # Define logscales
  for i = 1:length(logscale)
    ax[i], xlims[i] = set_log(plt[i], ax[i], xlims[i], logscale[i], "x",
                      :x, :xlerr, :xuerr, :set_xlim, :set_xscale)
    ax[i], ylims[i] = set_log(plt[i], ax[i], ylims[i], logscale[i], "y",
                      :y, :ylerr, :yuerr, :set_ylim, :set_yscale)
  end

  # Set axes limits
  ax1[:set_xlim](xlims[1]); ax1[:set_ylim](ylims[1])
  if ax_2  ax2[:set_xlim](xlims[2]); ax2[:set_ylim](ylims[2])  end

  # Set plot title
  ax1[:set_title](ti, fontsize=fntsiz+ti_offset)

  # Generate axes labels and legend, define axes label/tick colours
  ax1[:set_xlabel](xlab,fontsize=fntsiz+ax_offset)
  for n = 1:length(ylab)
    ax[n][:set_ylabel](ylab[n],fontsize=fntsiz+ax_offset, color=axcol[n])
    setp(ax[n][:get_yticklabels](),color=axcol[n])
  end

  if legpos[1] ≠ "None"
    ax1[:legend](fontsize=fntsiz+leg_offset, loc=legpos[1], ncol=legcol[1])
  end
  if ax_2 && legpos[2] ≠ "None"
    ax2[:legend](fontsize=fntsiz+leg_offset, loc=legpos[2], ncol=legcol[2])
  end

  # Set ticks and optional minor ticks
  if Mxtint > 0
    xint = collect(ax1[:get_xlim]()[1]:Mxtint:ax1[:get_xlim]()[2])
    ax1[:set_xticks](xint)
    if ax_2  ax2[:set_xticks](xint)  end
  end
  if Mytint[1] > 0  for i = 1:length(Mytint)
    yint = collect(ax[i][:get_ylim]()[1]:Mytint[i]:ax[i][:get_ylim]()[2])
    ax[i][:set_yticks](yint)
  end  end
  if mticks == "on"
    minorticks_on()
  else
    minorticks_off()
  end
  # Set minor x ticks
  if nmxt > 0
    mx = matplotlib[:ticker][:MultipleLocator](nmxt)
    for a in ax
      a[:xaxis][:set_minor_locator](mx)
    end
  end
  # Set minor y ticks
  for i = 1:length(nmyt)
    if nmyt[i] > 0
      my = matplotlib[:ticker][:MultipleLocator](nmyt[i])
      ax[i][:yaxis][:set_minor_locator](my)
    end
  end
  # Format ticks and frame
  Mtlen = tsc[1]⋅frw
  mtlen = tsc[2]⋅frw
  for a in ax
    a[:tick_params]("both", which="both", direction="in", top="on", right="on",
      labelsize=fntsiz, width=frw)
    a[:tick_params]("both", which="major", length=Mtlen)
    a[:tick_params]("both", which="minor", length=mtlen)
    a[:grid](linestyle=":", linewidth = frw)
    a[:spines]["bottom"][:set_linewidth](frw)
    a[:spines]["top"][:set_linewidth](frw)
    a[:spines]["left"][:set_linewidth](frw)
    a[:spines]["right"][:set_linewidth](frw)
  end
  tight_layout()

  # Return PyPlot data
  return fig, ax

end #function plot_data


###########################
###  PRIVATE FUNCTIONS  ###
###########################

### Functions associated with load_PlotData

"""
    rename_DF(pltdata, DFnames, err)

Rename columns in DataFrame `pltdata` with names specified in array `DFnames` with
entries for x, y, lower y error, upper y error, lower x error, and upper x error
based on the keyword given in `err`.

`err` (`String`):
- `"None"` (no errors, **default**)
- `"rangex"`, `"rangey"`, `"range"` (± x/y/x and y)
- `"pmrangex"`, `"pmrangey"`, `"pmrange"` (± x/y/x and y with different lower/upper values)
- `"percentx"`, `"percenty"`, `"percent"` (± err⋅x/y/x and y)
- `"pmpercentx"`, `"pmpercenty"`, `"pmpercent"` (± err⋅x/y/x and y with different lower/upper values)
- `"factorx"`, `"factory"`, `"factor"` (x/y/x and y ⋅1/err and ⋅err, respectively)
- `"pmfactorx"`, `"pmfactory"`, `"pmfactor"` (as above with different lower/upper values)
- `"valuex"`, `"valuey"`, `"value"` (err value directly taken from column)
"""
function rename_DF(pltdata, DFnames, err)
  if err == "None"
    names!(pltdata,DFnames[1:2])
  elseif startswith(err,"pm") && endswith(err,"x")
    names!(pltdata,DFnames[[1,2,5,6]])
  elseif endswith(err,"x")
    names!(pltdata,DFnames[[1,2,5]])
  elseif startswith(err,"pm") && endswith(err,"y")
    names!(pltdata,DFnames[1:4])
  elseif endswith(err,"y")
    names!(pltdata,DFnames[1:3])
  elseif startswith(err,"pm")
    names!(pltdata,DFnames)
  else
    names!(pltdata,DFnames[[1,2,3,5]])
  end

  return pltdata
end #function rename_DF


"""
    calc_errors(pltdata, col, err)

Calculate errors in DataFrame `pltdata` for columns with the specified column
names `col` based on the keyword `err`.

`err` (`String`):
- `"None"` (no errors, **default**)
- `"rangex"`, `"rangey"`, `"range"` (± x/y/x and y)
- `"pmrangex"`, `"pmrangey"`, `"pmrange"` (± x/y/x and y with different lower/upper values)
- `"percentx"`, `"percenty"`, `"percent"` (± err⋅x/y/x and y)
- `"pmpercentx"`, `"pmpercenty"`, `"pmpercent"` (± err⋅x/y/x and y with different lower/upper values)
- `"factorx"`, `"factory"`, `"factor"` (x/y/x and y ⋅1/err and ⋅err, respectively)
- `"pmfactorx"`, `"pmfactory"`, `"pmfactor"` (as above with different lower/upper values)
- `"valuex"`, `"valuey"`, `"value"` (err value directly taken from column)
"""
function calc_errors(pltdata, col, err)
  # ± errors (error range for x, y or both)
  # where `err` starts with `pm`, different values are
  # used for positive and negative errors
  if err == "rangex"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
  elseif err == "pmrangex"
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]]
  elseif err == "rangey"
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
  elseif err == "pmrangey"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]]
  elseif err == "range"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
  elseif err == "pmrange"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]]

  # ± (factor × value) (adding a percentage range)
  elseif err == "percentx"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
  elseif err == "pmpercentx"
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]].⋅pltdata[col[1]]
  elseif err == "percenty"
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
  elseif err == "pmpercenty"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]].⋅pltdata[col[2]]
  elseif err == "percent"
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
  elseif err == "pmpercent"
    pltdata[col[3]] = pltdata[col[2]] .- pltdata[col[3]].⋅pltdata[col[2]]
    pltdata[col[4]] = pltdata[col[2]] .+ pltdata[col[4]].⋅pltdata[col[2]]
    pltdata[col[5]] = pltdata[col[1]] .- pltdata[col[5]].⋅pltdata[col[1]]
    pltdata[col[6]] = pltdata[col[1]] .+ pltdata[col[6]].⋅pltdata[col[1]]

  # error × value (adding a factor error range)
  elseif err == "factorx"
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
  elseif err == "pmfactorx"
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[6]]
  elseif err == "factory"
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
  elseif err == "pmfactory"
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[4]]
  elseif err == "factor"
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[5]]
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[3]]
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
  elseif err == "pmfactor"
    pltdata[col[3]] = pltdata[col[2]] .⋅ 1.0./pltdata[col[3]]
    pltdata[col[4]] = pltdata[col[2]] .⋅ pltdata[col[4]]
    pltdata[col[5]] = pltdata[col[1]] .⋅ 1.0./pltdata[col[5]]
    pltdata[col[6]] = pltdata[col[1]] .⋅ pltdata[col[6]]
  end

  return pltdata
end #function calc_errors


### Functions associated with plot_data

"""
setup_axes(plot_list, twinax, ylab, logscale, xlims, ylims, Mytint, nmyt, cs, axcol, legpos, legcol)

If `twinax` is an array split `plot_list` into 2 datasets based on the indices
in `twinax`. Assure all paramters concerning y data are arrays of length 2, if
a second axis is allowed or length 1 otherwise. Return `plot_list` as array of
2 or 1 distinct lists an the adjusted parameters.
"""
function setup_axes(plot_list, twinax, ylab, logscale, xlims, ylims, Mytint, nmyt,
                    cs, axcol, legpos, legcol)

  # Set flag for second axis, asume no second axis
  ax_2 = false
  ax2 = nothing

  # Set up second axis, if the twinax array is defined
  plt1 = PlotData[]; plt2 = PlotData[]
  if !isempty(twinax)
    # Set flag true and define 2nd axis in PyPlot
    ax_2 = true
    ax2 = twinx()
    # Check correct input of twinax
    if length(twinax) ≠ length(plot_list)
      println("Array `twinax` must have the same length as array `phot_list`.")
      println("Script stopped."); exit()
    end

    # Assign data to the axes
    for i = 1:length(twinax)
      if twinax[i] == 1  push!(plt1,plot_list[i])
      else               push!(plt2,plot_list[i])
      end
    end

    # Make sure, all parameters for both axes are arrays of length 2
    if isa(logscale, String) logscale = String[logscale, logscale]  end
    if isa(axcol, String) axcol = String[axcol, axcol]  end
    if isa(cs, String) cs = String[cs, cs]  end
    if isa(Mytint, Int64) Mytint = Int64[Mytint, Mytint]  end
    if isa(nmyt, Int64) nmyt = Int64[nmyt, nmyt]  end
    if !isa(xlims, Array) xlims = [xlims, xlims]  end
    if !isa(ylims, Array) ylims = [ylims, ylims]  end
    if !isa(legpos, Array) legpos = [legpos, legpos]  end
    if isa(legcol, Int64) legcol = [legcol, legcol]  end
  else
    # Assign all data to axis 1, if no second axis
    plt1 = plot_list
    # If no 2nd axis, make sure, all parameters for both axes are arrays of length 1
    # and not single parameters (numbers, strings...)
    if isa(logscale, String) logscale = String[logscale]  end
    if isa(axcol, String) axcol = String[axcol]  end
    if isa(cs, String) cs = String[cs]  end
    if isa(Mytint, Int64) Mytint = Int64[Mytint]  end
    if isa(nmyt, Int64) nmyt = Int64[nmyt]  end
    if !isa(xlims, Array) xlims = [xlims]  end
    if !isa(ylims, Array) ylims = [ylims]  end
    if !isa(legpos, Array) legpos = [legpos]  end
    if isa(legcol, Int64) legcol = [legcol]  end
  end
  if isa(ylab, String)  ylab = String[ylab]  end
  if isempty(plt2)  plt = [plt1]
  else              plt = [plt1, plt2]
  end

  # Return adjusted data
  return plt, ax2, ax_2, ylab, logscale, xlims, ylims, Mytint, nmyt, cs, axcol, legpos, legcol
end #function setup_axes

function set_cs(plt, cs, lc, lt)
  for i = 1:length(cs)
    for j = 1:length(plt[i])
      if cs[i] == "own"
        plt[i][j].colour = lc[i][j]
        plt[i][j].dashes = lt[i][j]
      elseif cs[i] ≠ ""
        lc, dt = sel_ls(cs=cs[i], nc=j, nt=j)
        plt[i][j].colour = lc
        if plt[i][j].dashes[1] ≠ 0  plt[i].dashes = dt  end
      end
    end
  end

  return plt
end #function set_cs


"""
    plt_DataWithErrors(plt, ax, offset)

For each `PlotData` in array `plt` (and `ax`), retrieve the errors in the `PlotData` and
plot according to specifications. For error bars of markers, the standard cap size is `3`,
which can be adjusted by an `offset` (positive or negative number to be added to cap size).

Returns `fig` and `ax` (array of axes) for further plot modifications or printing.
"""
function plt_DataWithErrors(plt, ax, offset)
  # Redefine errors for error bars
  xerr = redef_err(plt,:x,:xlerr,:xuerr)
  yerr = redef_err(plt,:y,:ylerr,:yuerr)
  # Loop over graph data and plot each data according to its type and format
  # defined by the struct PlotData
  for i = 1:length(plt)
    if plt[i].yuerr ≠ nothing && plt[i].marker == "None"
      ax[:fill_between](plt[i].x, plt[i].ylerr, plt[i].yuerr,
          color=plt[i].colour, alpha=0.2)
      ax[:plot](plt[i].x, plt[i].y, lw = plt[i].lw,
          dashes=plt[i].dashes, color=plt[i].colour, label=plt[i].label)
    elseif plt[i].xuerr ≠ nothing && plt[i].yuerr ≠ nothing
      ax[:errorbar](plt[i].x, plt[i].y, xerr=[xerr[i][:lower], xerr[i][:upper]],
        yerr=[yerr[i][:lower], yerr[i][:upper]], lw = plt[i].lw,
        marker=plt[i].marker, dashes=plt[i].dashes, color=plt[i].colour,
        label=plt[i].label, capsize=3+offset)
    elseif plt[i].yuerr ≠ nothing
      ax[:errorbar](plt[i].x, plt[i].y, yerr=[yerr[i][:lower], yerr[i][:upper]],
        lw = plt[i].lw, marker=plt[i].marker, dashes=plt[i].dashes,
        color=plt[i].colour, label=plt[i].label, capsize=3+offset)
    elseif plt[i].xuerr ≠ nothing
      ax[:errorbar](plt[i].x, plt[i].y, xerr=[xerr[i][:lower], xerr[i][:upper]],
        lw = plt[i].lw, marker=plt[i].marker, dashes=plt[i].dashes,
        color=plt[i].colour, label=plt[i].label, capsize=3+offset)
    else
      ax[:plot](plt[i].x, plt[i].y,lw = plt[i].lw, marker=plt[i].marker,
      dashes=plt[i].dashes, color=plt[i].colour, label=plt[i].label)
    end
  end

  return plt, ax
end #function plt_DataWithErrors


"""
    redef_err(plt,val,low,high)

Recalculate errors relative to actual values (rather than absolute values in plot)
for use with error bars of markers. The function generates an array of DataFrames
with the revised errors in the columns `:upper` and `:lower`.
Data is generated from the array `plt` with the columns `val` with the measured/modelled
value and the columns with names `low` and `high` holding the
actual error values (rather than values relative to measured/modelled value.)
"""
function redef_err(plt,val,low,high)
  err = []
  for i = 1:length(plt)
    # Recalculate errors for error bars, if markers are plotted
    if plt[i].marker ≠ "None" && getfield(plt[i],high) ≠ nothing
      push!(err, DataFrame(upper = getfield(plt[i],high) .- getfield(plt[i],val),
            lower = getfield(plt[i],val) .- getfield(plt[i],low)))
    else
      # Otherwise, save errors as they are und column names `upper` and `lower`
      push!(err, DataFrame(upper = Float64[], lower = Float64[]))
    end
  end

  return err
end #function redef_err


"""
    set_log(plt, ax, lims, logscale, x, val, low, high, cmd1, cmd2)

Set x- and/or y-axis to log scale, if the string `logscale` consists of `"x"` or
`"y"` or `"xy"`. Adjust minimum/maximum values of the respective axis for logscale.
"""
function set_log(plt, ax, lims, logscale, x, val, low, high, cmd1, cmd2)
  if contains(logscale,x)
    if all([getfield(plt[i], high) == nothing for i =1:length(plt)])
      xmin = minimum(minimum.([getfield(plt[i], val) for i = 1:length(plt)]))
      xmax = maximum(maximum.([getfield(plt[i], val) for i = 1:length(plt)]))
    else
      xmin = minimum(minimum.([getfield(plt[i], low) for i = 1:length(plt)]))
      xmax = maximum(maximum.([getfield(plt[i], high) for i = 1:length(plt)]))
    end
    xmin = 10^floor(log10(xmin))
    xmax = 10^ceil(log10(xmax))
    lims = (xmin, xmax)
    ax[cmd1](lims)
    ax[cmd2]("log")
  end
  return ax, lims
end #function set_log


end #module pyp
