"""
# Module pltfcn

Routines to compile comparison plots TUV data and fitted functions
in a pdf in the output folder.

# Functions

- plot_j (public)
- beautify_chem (private)

"""
module pltfcn

#import modules and export functions
using DataFrames
using PyCall, PyPlot
export plot_j


"""
    plot_j(jvals, iofolder, fit, time)

Plot j values saved in dataframe jvals and parameterisations derived from parameters
in fit to file parameters.dat in the iofolder together with the time of creation (time).
"""
function plot_j(jvals,iofolder,fit,time)

  # Format titles with Latex using header of jvals
  ptitle = beautify_chem(names(jvals))

  # Initialise array of plots and define x data
  nj = length(fit)
  jplot = Array{PyCall.PyObject}(nj)
  xdata = jvals[1].*π./180.

  # Loop over all reactions
  for i=1:nj
    # get order of magnitude for plots
    p = floor(log10(fit[i].param[1]*exp(-fit[i].param[3])))
    # define parameters
    l = fit[i].param[1]/10^p
    m = fit[i].param[2]
    n = fit[i].param[3]

    # Clear plot and start new one
    clf()
    # Plot TUV data and fitted curve
    jplot[i] = scatter(xdata,jvals[i+1]./(10^p),color=:black,marker = "s",label="TUV data")
    jplot[i] = plot(xdata,l.*(cos.(xdata)).^m.*exp.(-n./cos.(xdata)), color=:red, label="Parameterisation")
    # Plot settings
    size(2400, 1800)
    xticks([0.,0.1*π,0.2*π,0.3*π,0.4*π,0.5*π],
           [L"$0$",L"$\frac{1}{10}\pi$",L"$\frac{1}{5}\pi$",L"$\frac{3}{10}\pi$",L"$\frac{2}{5}\pi$",L"$\frac{1}{2}\pi$"])
    # xticks(0:π/8:π/2)
    # ymax(5.e-5)
    xlabel("χ")
    xlab = latexstring("j\\ /\\ 10^{$(convert(Int64,p))}\\,\\mathrm{s}^{-1}")
    ylabel(xlab)
    xlim(0.0,π/2)
    ylim(0.0,ceil(jvals[1,i+1]/10^p))
    # axis(ymin=0.0)
    grid(linestyle="dotted")
    # tight_layout()
    ltitle = latexstring(ptitle[i+1])

    # Plot title, legend, and time of creation
    title(ltitle)
    legend()
    annotate("created $(Dates.format(time,"dd.mm.yyyy, HH:MM:SS"))",
    xy=[0.1;0.1],
    xycoords="axes fraction",
    ha="left",
    va="bottom")

    # save plot to temporary png
    savefig("$iofolder/j$i.png", dpi=300)
  end

  # compile pngs in single pdf and delete pngs
  compile_pdf(iofolder)
end #function plot_j

"""
    beautify_chem(reactions)

Format reactions with LaTeX for nicer titles in plots.
"""
function beautify_chem(reactions)

  # Initialise output
  chem = String[]
  # Loop over reactions and reformat with LaTeX
  for rxn in reactions
    lstr = convert(String,rxn) # convert dataframe symbols to strings
    lstr = replace(lstr, "->", "\\stackrel{h\\nu}{\\longrightarrow}") # format arrows
    lstr = replace(lstr, " + hv", "") # remove + hv (now above formatted arrows)
    lstr = replace(lstr, r"([A-Za-z)])([0-9]+)", s"\1_{\2}") # make numbers in chemical formulas subscripts: {\1}_{\2}
    lstr = replace(lstr, r"\(([0-9]+)", s"(^{\1}") # except for atom's energetic state
    lstr = replace(lstr, ".", "^{.}") # raise radical dots
    lstr = replace(lstr, r"([0-9]+) ", s"\1\\,") # insert halfspace between stoichiometric indices and species
    lstr = replace(lstr, r"=", "\\!=\\!") # no spaces between double bonds
    lstr = replace(lstr, r"([A-Za-z0-9])-([A-Za-z])", s"\1\\!-\\!\2") # format dashes
    lstr = replace(lstr, r"C:", "{\\ddot C}") # format biradical functions
    lstr = replace(lstr, r"CH:", "{\\ddot C}H") # format biradical functions
    # define and correct exceptions:
    lstr = replace(lstr, "further products", "\\text{further products}") # transform to text
    # Ensure unslanted font
    lstr = "\\mathrm{"*lstr*"}"

    # Save reformatted reaction to output
    push!(chem, lstr)
  end

  # Return reformatted output
  return chem
end #function beautify_chem


"""
    compile_pdf(iofolder)

Compile single png plots in one pdf file
and warn about possible missing software on failure.
"""
function compile_pdf(iofolder)

  try
    # compile pngs in a single pdf and delete all pngs
    run(`convert -quality 100 $iofolder/'*'.png $iofolder/plots.pdf`)
    cd(iofolder)
    pfiles = filter(r"./*\.png",readdir())
    for i = 1:length(pfiles)  rm(pfiles[i],force=true)  end
    cd("..")
  catch
    # if compilation fails issue warning about possible missing software
    println("\033[95mWARNING! UNIX tool 'convert' missing or unable to process plots.")
    println("No concatenation of plots in single pdf in results folder.")
    println("\033[96mMac users install imagemagick, e.g. using homebrew:")
    println("brew install imagemagick\033[0m")
  end
end #function compile_pdf
end #module fitfcn
