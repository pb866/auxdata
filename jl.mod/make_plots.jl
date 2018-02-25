__precompile__()


"""
# Module make_plots

Generate plots from DSMACC model output.

# Functions (public)
- lineplot
- plot_stack
- plot_flux
- plot_prodloss
- sel_ls
- night_shade
"""
module make_plots


##################
###  PREAMBLE  ###
##################

# Import external modules
try using PyPlot
catch
  Pkg.add("PyPlot")
  using PyPlot
end
try using DataFrames
catch
  Pkg.add("DataFrames")
  using DataFrames
end

# Export public functions
export lineplot,
       plot_stack,
       plot_flux,
       plot_prodloss,
       sel_ls,
       night_shade


###################
###  FUNCTIONS  ###
###################


"""
    lineplot(xdata,ydata,label,what,unit,icase,plotdata,nights,pltnight,t_frmt,sfile)

Generate PyPlot data for line plots from the model times specified in `xdata`
in the time format `t_frmt` and `ydata` with specified `label`s and in specified
`unit`s and return the raw plot data.

`what` specifies, whether to plot `"specs"` or `"rates"`, `icase` specifies the
current scenario(s), `plotdata` the species/reactions to be plotted.

Night-time periods with times specified in `nights` are plotted in the format
(color/transparancy) specified in `pltnight`.

If a file name `sfile` is provided, plots are saved to folder `FIG` as single pdfs
in addition to the pdf with the compiled plots as specified by the second script argument.

The function returns a PyCall.PyObject with the plot data.
"""
function lineplot(xdata,ydata,label,what,unit,icase,plotdata,nights,pltnight,t_frmt,sfile)

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
    if unit!="mlc"  &&  unit!="cm-3"
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
  # x axis
  if t_frmt == "TIME"
    xlim(xmin=xdata[1], xmax=xdata[end]+1.e-5)
    ax[:axes][:get_xaxis]()[:set_ticks](collect(0:12:xdata[end]+1.e-5))
    minorticks_on()
    mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
    ax[:xaxis][:set_minor_locator](mx)
    xlabel("model time / hours")
  elseif t_frmt == "JTIME"
    xlim(xmin=xdata[1], xmax=xdata[end])
    majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b, %H:%M")
    minorformatter = matplotlib[:dates][:DateFormatter]("")
    majorlocator = matplotlib[:dates][:HourLocator](byhour=(0, 12))
    minorlocator = matplotlib[:dates][:HourLocator](byhour=(3, 6, 9, 15, 18, 21))
    ax[:xaxis][:set_major_formatter](majorformatter)
    ax[:xaxis][:set_minor_formatter](minorformatter)
    ax[:xaxis][:set_major_locator](majorlocator)
    ax[:xaxis][:set_minor_locator](minorlocator)
    fig[:autofmt_xdate](bottom=0.2,rotation=-30,ha="left")
    xlabel("time (UTC)")
  end
  # y axis
  ymin, ymax = ylim()
  if ymin < 0  ylim(ymin=0)  end
  if unit=="mlc"  ||  unit=="cm-3"
    # Default y label
    if what=="specs" yl = "concentration / 10\$^{$(Int(p))}\$ molecules cm\$^{-3}\$"
    elseif what=="rates"  yl = "rate /cm\$^{3(n-1)}\$ molecules\$^{-n+1}\$ s\$^{-1}\$"
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

  # Night-time shading
  for n = 1:length(nights)÷2
    ax[:axvspan](xdata[nights[n,1]], xdata[nights[n,2]],
      facecolor=pltnight[1], alpha=pltnight[2])
  end

  # Save and return plot data
  if sfile != ""  savefig(sfile)  end
  close(fig)
  return fig
end #function lineplot


"""
    plot_stack(xdata,ylines,ystack,scenario,label,unit,lt,nights,pltnight,t_frmt,sfile)

Plot the concentrations (mainly of lumped sum species) in the specified `unit`
as a stacked area plots with boundary lines using the model time `xdata` in the
format `t_frmt` and `ylines` for the boundaries and `ystack` for the area in the
graphs.

The colour schemes and boundary line types are defined in tuple `lt` with the
colour scheme as String array as first entry and the line types with specifications
for `dashes` as second entry.

The `scenario` names are printed in the titel. Furthermore, `label`s are us for
the legend.

Night-time periods with times specified in `nights` are plotted in the format
(color/transparancy) specified in `pltnight`.

If a file name `sfile` is provided, plots are saved to folder `FIG` as single pdfs
in addition to the pdf with the compiled plots as specified by the second script argument.

The function returns a PyCall.PyObject with the plot data.
"""
function plot_stack(xdata,ylines,ystack,scenario,label,unit,lt,nights,pltnight,t_frmt,sfile)

  # Initialise current plot
  fig, ax = subplots()

  # Find y-axis maximum
  if unit!="mlc"  &&  unit!="cm-3"
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
  # x axis
  if t_frmt == "TIME"
    ax[:axes][:get_xaxis]()[:set_ticks](collect(0:12:xdata[end]+1.e-5))
    xlim(xmin=xdata[1], xmax=xdata[end]+1.e-5)
    minorticks_on()
    mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
    ax[:xaxis][:set_minor_locator](mx)
    xlabel("model time / hours")
  elseif t_frmt == "JTIME"
    xlim(xmin=xdata[1], xmax=xdata[end])
    majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b, %H:%M")
    minorformatter = matplotlib[:dates][:DateFormatter]("")
    majorlocator = matplotlib[:dates][:HourLocator](byhour=(0, 12))
    minorlocator = matplotlib[:dates][:HourLocator](byhour=(3, 6, 9, 15, 18, 21))
    ax[:xaxis][:set_major_formatter](majorformatter)
    ax[:xaxis][:set_minor_formatter](minorformatter)
    ax[:xaxis][:set_major_locator](majorlocator)
    ax[:xaxis][:set_minor_locator](minorlocator)
    fig[:autofmt_xdate](bottom=0.2,rotation=-30,ha="left")
    xlabel("time (UTC)")
  end

  # y axis
  ymin, ymax = ylim()
  if ymin < 0  ylim(ymin=0)  end
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

  # Night-time shading
  for n = 1:length(nights)÷2
    ax[:axvspan](xdata[nights[n,1]], xdata[nights[n,2]],
      facecolor=pltnight[1], alpha=pltnight[2])
  end

  # Save and return plot data
  if sfile != ""  savefig(sfile)  end
  close(fig)
  return fig
end #function plot_stack


"""
    plot_flux(spc, scenario, modtime, fluxes, unit, cs, nights, pltnight, t_frmt, sfile)

Plot `fluxes` (positive or negative) for species `spc` in the specified `scenario`
and the specified `unit` over model time (`modtime`) in the format `t_frmt` using
the colour scheme `cs` as time-resolved stacked area plot.

Night-time periods with times specified in `nights` are plotted in the format
(color/transparancy) specified in `pltnight`.

If a file name `sfile` is provided, plots are saved to folder `FIG` as single pdfs
in addition to the pdf with the compiled plots as specified by the second script argument.

The function returns a PyCall.PyObject with the plot data.
"""
function plot_flux(spc, scenario, modtime, fluxes, unit,
                   cs, nights, pltnight, t_frmt, sfile)
  # Generate stacked subplots for source and sink fluxes
  fig, ax = subplots()

  # Find maximum order of magnitude
  if unit!="mlc"  &&  unit!="cm-3"
    p = 0
  else
    p=floor(log10(maximum(abs.(sum(fluxes[1],1)))))
    if p==Inf || p==-Inf  p = 0  end
  end
  # Redesign data
  fluxes[1] .*= 10^-p

  # Fill plots with graphs
  ax[:stackplot](modtime,fluxes[1],colors=cs)
  # Define title and legend
  title("Scenario $scenario")
  ax[:legend](fluxes[2],fontsize=8)
  # Axes ticks, min/max, adjust maximum for rounding errors with +.00001
  # x axis
  if t_frmt == "TIME"
    ax[:axes][:get_xaxis]()[:set_ticks](collect(0:12:modtime[end]+1.e-5))
    xlim(xmin=modtime[1], xmax=modtime[end]+1.e-5)
    minorticks_on()
    mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
    ax[:xaxis][:set_minor_locator](mx)
    xlabel("model time / hours")
  elseif t_frmt == "JTIME"
    xlim(xmin=modtime[1], xmax=modtime[end])
    majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b, %H:%M")
    minorformatter = matplotlib[:dates][:DateFormatter]("")
    majorlocator = matplotlib[:dates][:HourLocator](byhour=(0, 12))
    minorlocator = matplotlib[:dates][:HourLocator](byhour=(3, 6, 9, 15, 18, 21))
    ax[:xaxis][:set_major_formatter](majorformatter)
    ax[:xaxis][:set_minor_formatter](minorformatter)
    ax[:xaxis][:set_major_locator](majorlocator)
    ax[:xaxis][:set_minor_locator](minorlocator)
    fig[:autofmt_xdate](bottom=0.2,rotation=-30,ha="left")
    xlabel("time (UTC)")
  end
  # y axis
  if unit=="mlc"  ||  unit=="cm-3"
    # Default y label
    ax[:set_ylabel]("$spc fluxes / \$10^{$(Int(p))}\$ molecules cm\$^{-3}\$ s\$^{-1}\$")
  else
    # Labels for converted units with volume mixing ratios
    ax[:set_ylabel]("$spc fluxes / $unit h\$^{-1}\$")
  end
  # annotate missing fluxes in the plot
  if fluxes[3] != "no fluxes"
    annotate("omitted fluxes:\n$(join(fluxes[3],'\n'))",
      xy=[2,0.95⋅ax[:get_ylim]()[2]], verticalalignment="top",fontsize=8)
  end
  # Define dotted grid lines
  ax[:grid](linestyle=":")
  # Night-time shading
  for n = 1:length(nights)÷2
    ax[:axvspan](modtime[nights[n,1]], modtime[nights[n,2]],
      facecolor=pltnight[1], alpha=pltnight[2])
  end
  # Final layout settings
  tight_layout()

  if sfile != ""  savefig(sfile)  end
  close(fig)
  # Reset data
  fluxes[1] ./= 10^-p
  # Return plot raw data for further pyplot processing
  return fig
end #function plot_flux


"""
    plot_prodloss(spc, scenario, modtime, source, sink,
                  unit, nights, pltnight, t_frmt, sfile)

Plot `source` (positive) and `sink` (negative) data for species `spc` in the
specified `scenario` and the specified `unit` over model time (`modtime`) in the
format `t_frmt` as time-resolved stacked area plot.

Night-time periods with times specified in `nights` are plotted in the format
(color/transparancy) specified in `pltnight`.

If a file name `sfile` is provided, plots are saved to folder `FIG` as single pdfs
in addition to the pdf with the compiled plots as specified by the second script
argument.

The function returns a PyCall.PyObject with the plot data.
"""
function plot_prodloss(spc, scenario, modtime, source, sink, unit,
                       nights, pltnight, t_frmt, sfile)
  # Generate stacked subplots for source and sink fluxes
  fig, (ax1, ax2) = subplots(2, sharex=true)
  # Define grid in both subplots
  ax1[:grid](linestyle=":"); ax2[:grid](linestyle=":")

  # Find maximum order of magnitude
  if unit!="mlc"  &&  unit!="cm-3"
    p = 0
  else
    p=floor(log10(min(maximum(sum(source[1],1)),maximum(abs.(sum(sink[1],1))))))
    if p==Inf || p==-Inf  p = 0  end
  end
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
  # x axis
  if t_frmt == "TIME"
    ax2[:axes][:get_xaxis]()[:set_ticks](collect(0:12:modtime[end]+1.e-5))
    xlim(xmin=modtime[1], xmax=modtime[end]+1.e-5)
    minorticks_on()
    mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
    ax2[:xaxis][:set_minor_locator](mx)
    xlabel("model time / hours")
  elseif t_frmt == "JTIME"
    xlim(xmin=modtime[1], xmax=modtime[end])
    majorformatter = matplotlib[:dates][:DateFormatter]("%d. %b, %H:%M")
    minorformatter = matplotlib[:dates][:DateFormatter]("")
    majorlocator = matplotlib[:dates][:HourLocator](byhour=(0, 12))
    minorlocator = matplotlib[:dates][:HourLocator](byhour=(3, 6, 9, 15, 18, 21))
    ax2[:xaxis][:set_major_formatter](majorformatter)
    ax2[:xaxis][:set_minor_formatter](minorformatter)
    ax2[:xaxis][:set_major_locator](majorlocator)
    ax2[:xaxis][:set_minor_locator](minorlocator)
    fig[:autofmt_xdate](bottom=0.2,rotation=-30,ha="left")
    xlabel("time (UTC)")
  end
  # y axis
  pextr=maximum(sum(source[1][:,2:end],1))
  nextr=abs(minimum(sum(sink[1][:,2:end],1)))
  extr = max(pextr,nextr)
  # p=10^floor(log10(extr))
  ymin = -ceil(extr)
  ymax = ceil(extr)
  ax1[:set_ylim](0,ymax); ax2[:set_ylim](ymin,0)
  if unit=="mlc"  ||  unit=="cm-3"
    # Default y label
    ax2[:set_ylabel]("$spc sink and source fluxes /\n\$10^{$(Int(p))}\$ molecules cm\$^{-3}\$ s\$^{-1}\$",ha="left")
  else
    # Labels for converted units with volume mixing ratios
    ax2[:set_ylabel]("$spc sink and source fluxes / $unit h\$^{-1}\$",ha="left")
  end
  # annotate missing fluxes in the plot
  if source[3] != "no fluxes"
    ax1[:annotate]("omitted sources:\n$(join(source[3],'\n'))",
      xy=[2,0.95⋅ax1[:get_ylim]()[2]], verticalalignment="top",fontsize=8)
  end
  if sink[3] != "no fluxes"
    ax2[:annotate]("ommitted sinks:\n$(join(sink[3],'\n'))",
      xy=[2,0.95⋅ax2[:get_ylim]()[1]],fontsize=8)
  end

  # Night-time shading
  for n = 1:length(nights)÷2
    ax1[:axvspan](modtime[nights[n,1]], modtime[nights[n,2]],
      facecolor=pltnight[1], alpha=pltnight[2])
    ax2[:axvspan](modtime[nights[n,1]], modtime[nights[n,2]],
      facecolor=pltnight[1], alpha=pltnight[2])
  end

  # Final layout settings
  tight_layout()
  subplots_adjust(hspace=0)

  if sfile != ""  savefig(sfile)  end
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
  src = ["#330099", "#9900CC", "#003DF5", "#7547FF", "#339D9E",
         "#33FFCC", "#66FF33", "#009933", "#998000", "#DBA708",
         "#99CC66", "#00FF80", "#00FFFF", "#6699CC", "#006699",
         "#0033FF", "#000088", "#274E13", "#A7B763", "#48BD14"]
  snk = ["#FF0000", "#FF8000", "#FFFF00", "#FF7A7A", "#FF0080",
         "#FF00FF", "#CC0099", "#CC0033", "#CC3300", "#FF470A",
         "#FFA347", "#FFCC99", "#FF9999", "#FF99CC", "#FFFF99",
         "#EEB702", "#AA8C2C", "#6D5A1C", "#B7AD8E", "#E69138"]
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


"""
    night_shade(jNO2)

Obtain an n×2 matrix with indices of the model times for the start/end of each night
from an array of `jNO2` values.
"""
function night_shade(jNO2)
  # Initilise arrays with indices for model times at sunset/sunrise
  sunset = Int64[]; sunrise = Int64[]
  # Loop over j(NO2) and find sunrises/sunsets by checking, whether j(NO2) changes to 0
  for i = 2:length(jNO2)-1
    if jNO2[i] == 0.0 && jNO2[i-1] > 0.0      push!(sunset,i)
    elseif jNO2[i] == 0.0 && jNO2[i+1] > 0.0  push!(sunrise,i)
    end
  end
  # If model run starts or ends at night add the respective indices for sunset/sunrise
  if all(jNO2[1:2] .== 0.0)  unshift!(sunset,1)  end
  if all(jNO2[end-1:end] .== 0.0)  push!(sunrise,length(jNO2))  end

  return return hcat(sunset,sunrise)
end #function night_shade

end #module make_plots
