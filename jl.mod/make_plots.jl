__precompile__()


"""
  # Module make_plots

  Generate plots from DSMACC model output
"""
module make_plots


  ##################
  ###  PREAMBLE  ###
  ##################

  # Import external modules
  using PyPlot, DataFrames

  # Export public functions
  export lineplot,
         compile_pdf,
         sel_ls


###################
###  FUNCTIONS  ###
###################


"""
    lineplot(xdata,ydata,label,what,unit,icase,plotdata,fname::String="DSMACC.pdf")

  Generate a line plots from `xdata` and `ydata` with specified `label`s
  and in specified `unit`s into the output file `fname`.

  `what` specifies, whether to plot `"specs"` or `"rates"`, `icase` specifies the
  current scenario(s), `plotdata` the species/reactions to be plotted.
"""
function lineplot(xdata,ydata,label,what,unit,
                  icase,plotdata,fname::String="DSMACC.pdf")#; unit::String="cm-3"

  # Initialise current plot
  fig, ax = subplots()
  maxval = [] #maximum value needed for axis formats
  ip = 0 #counter for y data needed to assign line style
  # Loop over different species/rates
  for n = 1:length(plotdata)
    # Loop over scenarios
    for m in icase
      # Retrieve maximum for current y dataset
      push!(maxval,maximum(ydata[m][Symbol(plotdata[n])]))
    end
  end
  # Retrieve order of overall maximum
  p=floor(log10(maximum(maxval)))
  factor = 10^-p

  # Loop over species/reactions and scenarios
  for spc in plotdata  for case in icase
    # Unit conversions
    if unit!="mlc"
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
    clr, lin = sel_ls(ip,ip) #select line styles
    # plot data
    plot(xdata,ydata[case][Symbol(spc)].*factor,label="$spc ($(label[case]))",
         color=clr, linestyle=lin)
  end  end

  # Format plots
  # Axes ticks, min/max
  ax[:axes][:get_xaxis]()[:set_ticks](collect(0:12:xdata[end]))
  ax[:set_xlim](xmin=xdata[1], xmax=xdata[end])
  minorticks_on()
  mx = matplotlib[:ticker][:MultipleLocator](3) # Define interval of minor ticks
  ax[:xaxis][:set_minor_locator](mx)
  ax[:set_ylim](ymin=0)
  # Axes labels
  xlabel("model time / hours")
  if unit=="mlc"
    # Default y label
    if what=="specs" yl = "concentration / 10\$^{$(Int(p))}\$ molecules cm\$^{-3}\$"
    elseif what=="rates"  yl = "rate /cm\$^{3\\cdot n}\$ molecules\$^{-n}\$ s\$^{-1}\$"
    end
  else
    if what=="specs" yl = "volume mixing ratio / $unit"
    elseif what=="rates"  yl = "rate / $unit\$^{-n+1}\$ h\$^{-1}\$"
    end
  end
  ylabel(yl)
  # Legend, grid, general layout
  legend() #ncol=2
  grid(linestyle=":")
  tight_layout()
  # Save plot to pdf
  savefig(fname)
end #function lineplot


"""
    compile_pdf(iofolder)

Compile single png plots in one pdf file
and warn about possible missing software on failure.
"""
function compile_pdf(iofolder::String=".",fname::String="plots",ext::String="pdf";
                     ndel::Int64=0)

  try
    # compile single plots in a one overall pdf and delete all single plots
    run(`convert -density 300 $iofolder/'*'.$ext -quality 100 $fname.pdf`)
    for i = 1:ndel  rm("$iofolder/$(lpad(i,4,"0")).$ext",force=true)  end
  catch
    # if compilation fails issue warning about possible missing software
    println("\033[95mWARNING! UNIX tool 'convert' missing or unable to process plots.")
    println("No concatenation of plots in single pdf in results folder.")
    println("\033[96mMac users install imagemagick, e.g. using homebrew:")
    println("brew install imagemagick\033[0m")
  end
end #function compile_pdf


"""
    sel_ls(nc::Int64=1,nt::Int64=1)

  Select line color and type from a preset list with up to 8 different line styles.
"""
function sel_ls(nc::Int64=1,nt::Int64=1)
  # Manually define line colours and styles
  lc = ["red","blue","green","#FF8533","#FF00FF","#9400D3","black","#A0A0A0"]
  lt = ["-","--","-.",":","-","--","-.",":"]

  # Return a set of colours/styles depending on the input choice
  return lc[nc], lt[nt]
end #function sel_ls

end #module JLplot
