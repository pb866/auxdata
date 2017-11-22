__precompile__()


"""
# Module fitfcn

Fit MCM parameterisations to TUV output:
    j / s^-1 = l·cos^m(χ)·exp(-n·sec(χ))
"""
module fitfcn

using DataFrames
using LsqFit
export fit_j


"""
    fit_j(jvals)

Derive MCM parameterisations for dataframe jvals with χ-dependent _j_ values.
    j / s^-1 = l·cos^m(χ)·exp(-n·sec(χ))
"""
function fit_j(jvals)

  # Define x values and transform from deg to rad
  xdata = jvals[:sza].data
  xdata = xdata.* π./180.
  # Initialise arrays for fitting data
  fit = []; sigma = []; rmse = []; R2 = []

  delrxn = Int64[]
  # Loop over all j data
  for i = 2:length(jvals)
    ydata = jvals[i].data # define y data
    if ydata[1] == 0.
      println("\033[95mReaction with no data skipped:\033[0m")
      println(convert(String,names(jvals)[i]))
      push!(delrxn,i)
      continue
    end

    # Define fit function with initial guesses
    model(x,p) = p[1].*(cos.(x)).^(p[2]).*exp.(-p[3]./cos.(x))
    p0 = [ydata[1],0.8,0.3]

    # Derive fit
    push!(fit, curve_fit(model, xdata, ydata, p0))
    # Derive sigma with 95% confidence
    push!(sigma, estimate_errors(fit[i-1-length(delrxn)],0.95))
    # Calculate statistical data for RMSE and R^2
    ss_err = sum(fit[i-1-length(delrxn)].resid.^2)
    ss_tot = sum((ydata-mean(ydata)).^2)
    # RMSE
    push!(rmse, √(ss_err/fit[i-1-length(delrxn)].dof))
    # R^2
    push!(R2, 1. - (ss_err/ss_tot))
  end
  for i in delrxn  delete!(jvals, i)  end

  return fit, sigma, rmse, R2

end #function fit_j

end #module fitfcn
