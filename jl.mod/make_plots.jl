__precompile__()


"""
# Module make_plots

Generate plots from DSMACC model output.

# Functions (public)
- lineplot
- plot_flux
- plot_prodloss
- sel_ls
"""
module make_plots


##################
###  PREAMBLE  ###
##################

# Import external modules
using PyPlot, DataFrames

# Export public functions
export lineplot,
       plot_stack,
       plot_flux,
       plot_prodloss,
       sel_ls


###################
###  FUNCTIONS  ###
###################


"""
    lineplot(xdata,ydata,label,what,unit,icase,plotdata)

Generate line plot data from `xdata` and `ydata` with specified `label`s
and in specified `unit`s and return the raw plot data.

`what` specifies, whether to plot `"specs"` or `"rates"`, `icase` specifies the
current scenario(s), `plotdata` the species/reactions to be plotted.

The function returns a PyCall.PyObject with the plot data.
"""
function lineplot(xdata,ydata,label,what,unit,icase,plotdata)

  # Initialise current plot
  fig, ax = subplots()
  maxval = [] #maximum value needed for axis formats
  ip = 0 #counter for y data needed to assign line style
  # Loop over different species/rates in all scenarios
  for n = 1:length(plotdata)  for m in icase
    # Retrieve maximum for current y dataset
    push!(maxval,maximum(ydata[m][Symbol(plotdata[n])]))
  end  end
  # Retrieve order of overall maximum
  p=floor(log10(maximum(maxval)))
  factor = 10^-p

  # Loop over species/reactions and scenarios
  for spc in plotdata  for case in icase
    # Unit conversions
    if unit!="mlc"  ||  unit!="cm-3"
      p = 0
      factor = 1./ydata[case][:M]
      if what == "rates"  factor *= 3600.  end
    end
    if unit=="ppm"
      factor *= 1.e6
    elseif unit=="ppb"
      factor *= 1.e9
    elseif unit=="ppt"
      factor *= 1.e12
    end
    # Generate plots
    ip += 1
    clr, dt = sel_ls(cs="line",nc=ip,nt=ip) #select line styles
    # plot data
    plot(xdata,ydata[case][Symbol(spc)].*factor,label="$spc ($(label[case]))",
         color=clr, dashes=dt)
  end  end

  # Format plots
  # Axes ticks, min/max
  ax[:axes][:get_xaxis]()[:set_ticks](collect(0:12:xdata[end]))
  xlim(xmin=xdata[1], xmax=xdata[end])
  minorticks_on()
  mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
  ax[:xaxis][:set_minor_locator](mx)
  ymin, ymax = ylim()
  if ymin < 0  ylim(ymin=0)  end
  # Axes labels
  xlabel("model time / hours")
  if unit=="mlc"  ||  unit=="cm-3"
    # Default y label
    if what=="specs" yl = "concentration / 10\$^{$(Int(p))}\$ molecules cm\$^{-3}\$"
    elseif what=="rates"  yl = "rate /cm\$^{3\\cdot n}\$ molecules\$^{-n}\$ s\$^{-1}\$"
    end
  else
    # Labels for converted units with volume mixing ratios
    if what=="specs" yl = "volume mixing ratio / $unit"
    elseif what=="rates"  yl = "rate / $unit\$^{-n+1}\$ h\$^{-1}\$"
    end
  end
  ylabel(yl)
  # Legend, grid, general layout
  legend() #ncol=2
  grid(linestyle=":")
  tight_layout()
  # Save and return plot data
  close(fig)
  return fig
end #function lineplot


"""
    plot_stack(xdata,ylines,ystack,scenario,label,unit,lt)

Plot the concentrations (mainly of lumped sum species) as a stacked area plots
with boundaries using `xdata` and `ylines` for the boundaries and `ystack` for
the area in the graphs.

The colour schemes and boundary line types are defined in tuple `lt` with the
colour scheme as String array as first entry and the line types with specifications
for `dashes` as second entry.

The `scenario` names are printed in the titel. Furthermore, `label`s are us for
the legend and `unit` is used for the y axis label.

The function returns a PyCall.PyObject with the plot data.
"""
function plot_stack(xdata,ylines,ystack,scenario,label,unit,lt)

  # Initialise current plot
  fig, ax = subplots()

  # Find y-axis maximum
  if unit!="mlc"  ||  unit!="cm-3"
    p = 0
  else
    p=floor(log10(maximum(ylines[end])))
    if p==Inf || p==-Inf  p = 0  end
  end
  # Redesign data
  ystack .*= 10^-p; ylines .*= 10^-p

  # Fill plots with graphs
  ax[:stackplot](xdata,ystack,colors=lt[1],alpha=0.3)
  ax[:legend](label,fontsize=8) #ncol=2
  for l = 1:length(ylines)
    ax[:plot](xdata,ylines[l],color=lt[1][l],dashes=lt[2])
  end

  ### Format plots ###
  # Axes ticks, min/max
  ax[:axes][:get_xaxis]()[:set_ticks](collect(0:12:xdata[end]))
  xlim(xmin=xdata[1], xmax=xdata[end])
  minorticks_on()
  mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
  ax[:xaxis][:set_minor_locator](mx)
  ymin, ymax = ylim()
  if ymin < 0  ylim(ymin=0)  end
  # Axes labels
  xlabel("model time / hours")
  if unit=="mlc"  ||  unit=="cm-3"
    # Default y label
    yl = "concentration / 10\$^{$(Int(p))}\$ molecules cm\$^{-3}\$"
  else
    # Labels for converted units with volume mixing ratios
    yl = "volume mixing ratio / $unit"
  end
  ylabel(yl)
  # Legend, grid, general layout
  title("Scenario $scenario")
  ax[:grid](linestyle=":")
  tight_layout()
  # Save and return plot data
  close(fig)
  return fig

end #function plot_stack


"""
    plot_flux(spc, scenario, modtime, fluxes, cs)

Plot `fluxes` (positive or negative) over model time (`modtime`) using the colour
scheme `cs` as time-resolved stacked area plot labelling the current species `spc`
in the y axis and the current `scenario` in the title.

The function returns a PyCall.PyObject with the plot data.
"""
function plot_flux(spc, scenario, modtime, fluxes, cs)
  # Generate stacked subplots for source and sink fluxes
  fig, ax = subplots()
  # Define dotted grid lines
  ax[:grid](linestyle=":")

  # Find maximum
  p=floor(log10(maximum(sum(fluxes[1],1))))
  if p==Inf || p==-Inf  p = 0  end
  # Redesign data
  fluxes[1] .*= 10^-p


  # Fill plots with graphs
  ax[:stackplot](modtime,fluxes[1],colors=cs)
  # Define title and legend
  title("Scenario $scenario")
  ax[:legend](fluxes[2],fontsize=8)
  # Axes ticks, min/max, adjust maximum for rounding errors with +.00001
  ax[:axes][:get_xaxis]()[:set_ticks](collect(0:12:modtime[end]+.00001))
  xlim(xmin=modtime[1], xmax=modtime[end]+.00001)
  minorticks_on()
  mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
  ax[:xaxis][:set_minor_locator](mx)
  # Axes labels
  xlabel("model time / hours")
  ax[:set_ylabel]("$spc fluxes / \$10^{$(Int(p))}\$ molecules cm\$^{-3}\$ s\$^{-1}\$")
  # annotate missing fluxes in the plot
  if fluxes[3] != "no fluxes"
    annotate("omitted fluxes:\n$(join(fluxes[3],'\n'))",
      xy=[2,0.95⋅ax[:get_ylim]()[2]], verticalalignment="top",fontsize=8)
  end
  # Final layout settings
  tight_layout()

  close(fig)
  # Reset data
  fluxes[1] ./= 10^-p
  # Return plot raw data for further pyplot processing
  return fig
end #function plot_flux


"""
    plot_prodloss(spc, scenario, modtime, source, sink)

Plot `source` (positive) and `sink` (negative) data over model time (`modtime`)
as time-resolved stacked area plot labelling the current species `spc`
in the y axis and the current `scenario` in the title.

The function returns a PyCall.PyObject with the plot data.
"""
function plot_prodloss(spc, scenario, modtime, source, sink)
  # Generate stacked subplots for source and sink fluxes
  fig, (ax1, ax2) = subplots(2, sharex=true)
  # Define grid in both subplots
  ax1[:grid](linestyle=":"); ax2[:grid](linestyle=":")

  # Find maximum
  p=floor(log10(min(maximum(sum(source[1],1)),maximum(sum(sink[1],1)))))
  # Redesign data
  source[1] .*= 10^-p; sink[1] .*= -10^-p


  # Set colour schemes
  cs_src, dt = sel_ls(cs="source",nc=1:length(source[2]))
  cs_snk, dt = sel_ls(cs="sink",nc=1:length(sink[2]))
  # Fill plots with graphs
  ax1[:stackplot](modtime,source[1],colors=cs_src)
  ax2[:stackplot](modtime,sink[1],colors=cs_snk)

  # Define title and legend
  ax1[:set_title]("Scenario $scenario")
  ax1[:legend](source[2],fontsize=8); ax2[:legend](sink[2],fontsize=8)
  # Axes ticks, min/max,\n adjust maximum for rounding errors with +.00001
  ax2[:axes][:get_xaxis]()[:set_ticks](collect(0:12:modtime[end]+.00001))
  xlim(xmin=modtime[1], xmax=modtime[end]+.00001)
  pextr=maximum(sum(source[1][:,2:end],1))
  nextr=abs(minimum(sum(sink[1][:,2:end],1)))
  pp=10^floor(log10(pextr))
  pn=10^floor(log10(nextr))
  ymin = -ceil(nextr/pn)⋅pn
  ymax = ceil(pextr/pp)⋅pp
  minorticks_on()
  mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
  ax2[:xaxis][:set_minor_locator](mx)
  ax1[:set_ylim](0,ymax); ax2[:set_ylim](ymin,0)
  # Axes labels
  xlabel("model time / hours")
  ax2[:set_ylabel]("$spc sink and source fluxes /\n\$10^{$(Int(p))}\$ molecules cm\$^{-3}\$ s\$^{-1}\$",ha="left")
  # annotate missing fluxes in the plot
  if source[3] != "no fluxes"
    ax1[:annotate]("omitted sources:\n$(join(source[3],'\n'))",
      xy=[2,0.95⋅ax1[:get_ylim]()[2]], verticalalignment="top",fontsize=8)
  end
  if sink[3] != "no fluxes"
    ax2[:annotate]("ommitted sinks:\n$(join(sink[3],'\n'))",
      xy=[2,0.95⋅ax2[:get_ylim]()[1]],fontsize=8)
  end
  # Final layout settings
  tight_layout()
  subplots_adjust(hspace=0)

  close(fig)
  # Reset data
  source[1] ./= 10^-p; sink[1] ./= -10^-p
  # Return plot raw data for further pyplot processing
  return fig
end #function plot_flux


"""
    sel_ls(;cs::String="line",nc=1,nt=1)

Select colours with index `nc` from a colour scheme with the key word `cs`
(currently `"line"`, `"sink"`, and `"source"` available) and a line type
(from a range of solid, dashed, dash-dotted, and dotted lines) with the index `nt`.

Single colours/line types using integers or subsets using unit ranges `m:n` may be used.
"""
function sel_ls(;cs::String="line",nc=1,nt=1)
  # Manually define sets of colour schemes
  lin = ["black", "red", "blue", "green", "#FF8533",
         "#D499ED", "#AC6012", "#FF00FF", "#A50029", "#808000",
         "#A0A0A0", "#FF7F7F", "#3D9EFF", "#99CC99", "#FBC96D",
         "#7F7FFF", "#DDBFA0", "#FFB2FF", "#C68696", "#D8D8B2"]
  src = ["#9900CC", "#003DF5", "#7547FF", "#33FFCC", "#66FF33", "#009933", "#998000",
         "#990066", "#330099", "#006699", "#00FFFF", "#00FF80", "#99CC66", "#6699CC",
         "#9966CC", "#CC66CC", "#B83DB8", "#8A2E8A", "#0033FF", "#000088"]
  snk = ["#FF0000", "#FF8000", "#FFFF00", "#FF7A7A", "#FF0080", "#FF00FF", "#CC0099",
         "#CC0033", "#CC3300", "#FF470A", "#FFA347", "#FFCC99", "#FF9999", "#FF99CC",
         "#FFFF99", "#EEB702", "#AA8C2C", "#6D5A1C", "#B7AD8E", "#E69138"]
  # Combine colour schemes in a dataframe
  colourschemes = DataFrame(line=lin, source=src, sink=snk)
  # Save selected colour or subset of colours
  lc = colourschemes[Symbol(cs)][nc]
  # Manually define line types
  dt = [[],[3,2,3,2],[4,1.6,1,1.6,1,1.6],[6,2,1,1,1,1,1,2],[1,1.6,1,1.6],
        [3,1.6,3,1.6,1,1.6],[4,1.6,1,1.6,1,1.6,4,1.6,1,1.6,1,1.6,1,1.6],
        [4,1.6,4,1.6,1,1.6,1,1.6,4,1.6,4,1.6,1,1.6],
        [4,1.6,4,1.6,4,1.6,1,1.6,1,1.6],[5,1.6,2,1.6],
        [],[3,2,3,2],[4,1.6,1,1.6,1,1.6],[6,2,1,1,1,1,1,2],[1,1.6,1,1.6],
        [3,1.6,3,1.6,1,1.6],[4,1.6,1,1.6,1,1.6,4,1.6,1,1.6,1,1.6,1,1.6],
        [4,1.6,4,1.6,1,1.6,1,1.6,4,1.6,4,1.6,1,1.6],
        [4,1.6,4,1.6,4,1.6,1,1.6,1,1.6],[5,1.6,2,1.6]]

  # Return a set of colours/styles depending on the input choice
  return lc, dt[nt]
end #function sel_ls

end #module make_plots
