__precompile__()

"""
# Module pyp

Load data and plot with PyPlot.


# Functions

- rd_data
- lineplot
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

# Load further modules
if all(LOAD_PATH.!=Base.source_dir())  push!(LOAD_PATH,Base.source_dir())  end
using fhandle: test_file
using make_plots: sel_ls
# Export public functions
export rd_data,
       lineplot


###################
###  FUNCTIONS  ###
###################

"""
    rd_data(ifile; ix::Int64=1, iy=0, headers=false, SF=1, sep::String="")

Read data from text file `ifile` in the following format:
    # Use '#' as first character for comments or data to be ignored.
    # Optional: specify column header names for DataFrame keywords
    # with keyword `jlheaders` as:
    # jlheaders: x1, y1, y2, ..., yn
    # (You may use whitespace, commas (,), semicolons (;), or pipes (|) as
    #  separators in the list above)
    <data separated by whitespace or any charecter/string>

The function uses several keyword arguments (**kwargs) for more freedom in the
file format or the selection of data.

### **kwargs

- `ix`: Column index for column in `ifile` holding the x data (default index: `1`;
  default column name in output DataFrame: `x`)
- `iy`: Column index/indices of y data columns in `ifile` (default index: `1`;
  default column name(s) in output DataFrame: `y1` ... `yn`)
  Columns can be specified using an integer for the selection of a single column,
  ranges (`<n(low)>:<n(high)>`), or arrays with special selections (`[n1, n2, ..., nn]`)
  where column order can be rearranged
- `headers`: If headers is set to `true` (default: `false`), you need to specify
  column header names for all columns for the output DataFrame as described above
  using the keyword `jlheaders`. All columns will be read in and saved the same
  order as in `ifile`.
- `SF`: You can optionally apply scaling factors to y data (default value: `1`
  (no scaling)). If `SF` is an integer, the scaling factor will be applied to all
  y columns. You can apply scaling to each y column individually by providing an
  array of scaling factors of length number of columns - 1 (for the x column). If
  you only want to scale certain column(s), set the scaling factors for these columns
  and use `1` otherwise.
- `sep`: You can specify any column separator with the keyword charactar `sep`
  (default: whitespace). Separators can be any unicode character (even special characters
  such as `≠` or `α`) or string series of unicode characters (including whitespace).
"""
function rd_data(ifile; ix::Int64=1, iy=0, headers=false, SF=1, sep::String="")
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
      ix = 1; iy = 0
    end
    # Find and delete comment lines
    del=find(startswith.(lines,"#"))
    deleteat!(lines,del)
  end

  # Determine number of y columns for default case
  if iy == 0
    iy = 2:length(split(lines[1]))
  end
  # Initilise x and y data
  x = Float64[]; y = Matrix{Float64}(0, length(iy))
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
    # Save current line to respective data arrays
    push!(x,float(raw[ix]))
    y = vcat(y,transpose(float.(raw[iy]).⋅SF))
  end

  # Generate output DataFrame
  output = DataFrame(x = x)
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
    lineplot(pdata, *args; **kwargs)
Plot `pdata` from a DataFrame to a pdf using PyPlot with the specifications
given by `*args` and `**kwargs`.


# *args

- `xlab`: String for x axis label (default: `model time / hours`)
- `ylab`: String for y axis label (default: `concentration / mlc cm\$^{-3}\$ s\$^{-1}\$`)


# **kwargs

- `cs`: Set colour scheme from function sel_ls, 3 schemes available:
  `"line"` (default), `"source"`, and `"sink"`
- `lw`: Set line width of graphs
- `lt`: Set line style of graphs (different set of solid, dashed, dash-dotted,
  and dotted styles available)
- `mticks`: Use string `"on"` to show minor ticks and `"off"` suppress display
  of minor ticks
- `nmxt`: Set interval value of minor x ticks
- `nmyt`: Set interval value of minor y ticks
- `tsc`: Tuple with scaling factors for length of major (first tuple entry) and
  minor (second tuple entry) ticks in relation to their width
- `xlims`: Tuple with minimum and maximum x values
- `ylims`: Tuple with minimum and maximum y values
- `figsiz`: Tuple with (x, y) dimensions in inch for figure size
- `fntsiz`: Value of font size for tick labels and legend, axis labels are increased
  by 2
- `frw`: Set line width of plot border lines, grid lines, and tick width.
"""
function lineplot(pdata, xlab::String="model time / hours",
                  ylab::String="concentration / mlc cm\$^{-3}\$ s\$^{-1}\$";
                  cs="line", lw::Number=1.4, lt=nothing, mticks::String="on",
                  nmxt::Number=0, nmyt::Number=0, xlims=nothing, ylims=nothing,
                  figsiz::Tuple{Number,Number}=(6,4), fntsiz::Number=12,
                  frw::Number=1, tsc::Tuple{Number,Number}=(4.5,2.5))
  # Start plot
  fig, ax = subplots(figsize=figsiz)
  # Add empty label column to input matrix, if missing
  if length(pdata[1,:]) == 2
    lab = String[]
    for i = 1:length(pdata[1,:])  push!(lab,"")  end
    pdata = hcat(pdata,lab)
  end

  # Set line colours and types
  nc = 1:length(pdata[:,1])
  if lt == nothing  lt = 1:length(pdata[:,1])  end
  lstyle = sel_ls(cs = cs, nc = nc, nt = lt)

  # Plot graphs
  for i = 1:length(pdata)÷3
    ax[:plot](pdata[i,1], pdata[i,2], linewidth=lw, dashes=lstyle[2][i],
              color=lstyle[1][i], label=pdata[i,3])
  end

  # Set boundaries
  xlim(xlims); ylim(ylims);

  # Generate axes labels and legend
  xlabel(xlab,fontsize=fntsiz+2); ylabel(ylab,fontsize=fntsiz+2)
  ax[:legend](fontsize=fntsiz)

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

  # close(fig)

  return fig
end #function lineplot

end #module pyp
