__precompile__()


"""
# Module VectorMaths

Data Analysis and Manipulation tools containing functions:

- `DFvector_mean` to return a DataFrame with common x data, adjusted y data,
  the mean and sum of the y data of two DataFrames
- `multivector_mean` to return a DataFrame with the common x data, adjusted y data,
  the mean and sum of the y data of any number of individaul DataFrames with an
  x and a y column either for their common x data range or over the whole x data
  range
"""
module VectorMaths


##################
###  PREAMBLE  ###
##################

# Load packages
try using DataFrames
catch
  Pkg.add("DataFrames")
  using DataFrames
end
try using Dierckx
catch
  Pkg.add("Dierckx")
  using Dierckx
end

export DFvector_mean,
       multivector_mean


##########################
###  PUBLIC FUNCTIONS  ###
##########################


"""
    DFvector_mean(df1::DataFrames.DataFrame, df2::DataFrames.DataFrame;
      x1::Symbol, y1::Symbol, x2::Symbol, y2::Symbol)

Calculate the mean of `y1(x1)` and `y2(x2)` from 2 different DataFrames `df1`
and `df2`. If `x1` and `x2` differ, combine both vectors in the overlapping range.

For missing y values in either `df1` or `df2`, cubic splines are used to interpolate
the data and have 2 consistent datasets with shared x values.

The function uses keyword arguments `x1`, `x2`, `y1`, and `y2` to specify the column
symbols in each DataFrame. If the symbols are obsolete, column 1 is assumed to hold
the x data and column 2 the y data in each DataFrame, respectively.

The function returns a DataFrames with columns `x` with the combined x values from
both datasets, columns `y1` and `y2` with the adjusted y data using the combined x
data of equal array length, a column `mean` with the mean of `y1` and `y2`, and
a column `sum` with the sum of `y1` and `y2`.
"""
function DFvector_mean(df1::DataFrames.DataFrame, df2::DataFrames.DataFrame;
  x1::Symbol=Symbol(""), y1::Symbol=Symbol(""), x2::Symbol=Symbol(""), y2::Symbol=Symbol(""))
  # Set default columns for x/y data
  if x1 == Symbol("")  x1 = names(df1)[1]  end
  if y1 == Symbol("")  y1 = names(df1)[2]  end
  if x2 == Symbol("")  x2 = names(df2)[1]  end
  if y2 == Symbol("")  y2 = names(df2)[2]  end
  # Find overlapping range in x data
  xmin = max(df1[x1][1],df2[x2][1])
  xmax = min(df1[x1][end],df2[x2][end])

  # Abort for non-strictly ascending x data
  if unique(sort(df1[x1])) ≠ df1[x1] || unique(sort(df2[x2])) ≠ df2[x2]
    println("No strictly ascending x data. Script stopped."); exit()
  end

  # Create cubic spline of both datasets incase to interpolate missing data
  spl1 = Spline1D(df1[x1], df1[y1])
  spl2 = Spline1D(df2[x2], df2[y2])

  # Combine x data of both DataFrames
  xdata = unique(sort(vcat(df1[x1],df2[x2])))
  imin = find(xdata.==xmin)[1]; imax = find(xdata.==xmax)[1]
  xdata = xdata[imin:imax]

  # Assing y data from each DataFrame, use cubic spines for missing data
  ydat1 = Number[]; ydat2 = Number[]
  for x in xdata
    try push!(ydat1,df1[y1][find(df1[x1].==x)[1]])
    catch; push!(ydat1,spl1(x))
    end
    try push!(ydat2,df2[y2][find(df2[x1].==x)[1]])
    catch; push!(ydat2,spl2(x))
    end
  end

  # Calculate sum and mean of both y-vectors
  ysum = (ydat1.+ydat2)
  ymean = ysum./2
  # Return combined x data, mean of combined y data,
  # and adjusted y data vectors of each vector in a DataFrame
  return DataFrame(x = xdata, y1 = ydat1, y2 = ydat2, mean = ymean, sum = ysum)
end #function DFvector_mean


"""
    multivector_mean(vectors::Vector{DataFrames.DataFrame}; datarange::String="all")

From a vector of DataFrames (`vectors`) with x values in the first column and
y values in the second column, return a single DataFrame with unified x data in
the first column `x`, the ydata in the following columns `y1` to `yn`, the mean
in column `mean` and the sum in column `sum`. The x data in each DataFrame must
be strictly monotonic.

There are two options for the keyword argument `datarange`: `common` and `all`.
If `common` is chosen, data is only returned for the overlapping range of all
vectors. If the data have different step sizes, all x values are collected and
cubic splines are used to interpolate missing data.

If `all` is chosen, all vectors are filled with `NaN`s, if their range is shorter
than the maximum range of all vectors. `NaN`s are ignored for the calculation of
the `sum` and the `mean`. Furthermore, for vectors of different lengths, the mean
is multiplied after each premature end of a vector by a scaling factor to assure
a contineous smooth line even if an outlier has ended:

    SF = Ø(y values of all vectors)/Ø(y values of ending vectors)

This means, that the edges of each shorter vector will influence the mean outside
the range of this vector.
"""
function multivector_mean(vectors::Vector{DataFrames.DataFrame}; datarange::String="all")

  ### Check input data
  test_monotonicity(vectors)

  ### Find common range of all vectors
  # Find start/end x value of each vector (DataFrame)
  mins = Float64[]; maxs = Float64[]
  for d in vectors
    push!(mins,d[1][1])
    push!(maxs,d[1][end])
  end
  # Find boundaries of common range in all datasets
  lowest_common  = maximum(mins)
  largest_common = minimum(maxs)

  ### Get all data points in each DataFrame within the common range
  common_xdata = get_xdata(vectors, lowest_common, largest_common)

  # Get cubic splines of all vectors
  spl = Dierckx.Spline1D[]
  for d in vectors  push!(spl,Spline1D(d[1],d[2]))  end

  # Assign y data from each DataFrame, use cubic spines for missing data
  common_ydata = get_ydata(common_xdata,vectors,mins,maxs,spl)


  # Generate output DataFrame with common x data, individual y data columns,
  # the mean and sum
  if datarange == "common"
    # Compile output DataFrame
    dfr = DataFrame(x = common_xdata)
    for i = 1:length(vectors)
      dfr[Symbol("y$i")] = common_ydata[:,i]
    end
    dfr[:mean] = vec(mean(common_ydata,2))
    dfr[:sum]  = vec(sum(common_ydata,2))
  else
    # Find all x and y data outside common range
    # Assign NaN's to missing data
    lower_xdata = get_xdata(vectors, 0, lowest_common, bounds="exclude")
    upper_xdata = get_xdata(vectors, largest_common, maximum(maxs), bounds="exclude")

    # Find indices in arrays, where some data columns end in the lower or upper data
    mi, Mi = find_boundaries(mins,maxs,lower_xdata,upper_xdata)

    lower_ydata = get_ydata(lower_xdata,vectors,mins,maxs,spl)
    upper_ydata = get_ydata(upper_xdata,vectors,mins,maxs,spl)

    SF = get_ScalingFactor(lower_xdata, common_xdata, upper_xdata,
                           lower_ydata, common_ydata, upper_ydata, mi, Mi)

    # Get mean and sum
    lower_mean, lower_sum = calc_MeanSum(lower_ydata,SF[1])
    upper_mean, upper_sum = calc_MeanSum(upper_ydata,SF[2])

    # Compile output DataFrame
    dfr = DataFrame(x = vcat(lower_xdata, common_xdata, upper_xdata))
    for i = 1:length(vectors)
      dfr[Symbol("y$i")] = vcat(lower_ydata[:,i], common_ydata[:,i], upper_ydata[:,i])
    end
    dfr[:mean] = vcat(lower_mean, vec(mean(common_ydata,2)), upper_mean)
    dfr[:sum]  = vcat(lower_sum, vec(sum(common_ydata,2)), upper_sum)
  end


  return dfr
end #function multivector_mean


###########################
###  PRIVATE FUNCTIONS  ###
###########################


"""
    test_monotonicity(vectors)

Test the the first column in each DataFrame in a vector of DataFrames (`vectors`)
is strictly monotonic. Issue a warning and stop the script, if not.
"""
function test_monotonicity(vectors)
  for v = 1:length(vectors)
    if vectors[v][1] != sort(unique(vectors[v][1]))
      println("Error! X data in dataframe $v not strictly monotonic.\n Script stopped.")
      exit(99)
    end
  end
end #function test_monotonicity

"""
    get_xdata(vectors,lb,ub;bounds::String="include")

For a vector of DataFrames (`vectors`), extract the x data from each first column
within the given bounds `lb` and `ub`. If the keyword argument `bounds` is set to
`"include"`, bounds are included, if set to `"exclude"`, bounds are excluded.
"""
function get_xdata(vectors,lb,ub;bounds::String="include")
  # Retrieve indices of all x datapoints within the specified boundaries
  vrange = []
  for d in vectors
    if bounds == "include"
      push!(vrange,find(lb.≤d[1].≤ub))
    else bounds == "exclude"
      push!(vrange,find(lb.<d[1].<ub))
    end
  end

  # Collect all possible datapoints from each vector and unify the data
  xdata = Float64[]
  for i = 1:length(vectors)
    xdata = vcat(xdata,vectors[i][1][vrange[i]])
  end
  xdata = unique(sort(xdata))

  # Return a unified x data array
  return xdata
end #function get_xdata


"""
    get_ydata(xdata,vectors,mins,maxs,spl)

From the unified `xdata` in a specified range in `vectors` (vector holding
DataFrames with x data in the first column and y data in the second column),
the minima (`mins`) and maxima (`maxs`) of each individual x data and an array
of cubic splines of each DataFrame (`spl`), return a matrix with unified ydata,
where missing datapoints due to different step sizes within the range of each
DataFrame are interpolated by a cubic spline in `spl` and are filled with `NaN`s
outside the range of the individual datasets.
"""
function get_ydata(xdata,vectors,mins,maxs,spl)
  # Initialise output matrix
  ydata = Matrix{Float64}(length(xdata),0)
  # Loop over all DataFrames
  for i = 1:length(vectors)
    # Initialise y data of current DataFrame
    yd = Float64[]
    # Loop over unified x data
    for x in xdata
      # Try to add value for current x value, if y is missing
      # interpolate with cubic spline, if inside the range of the current vector
      # or fill with NaN outside its range
      try push!(yd,vectors[i][2][find(vectors[i][1].==x)[1]])
      catch
        if mins[i] ≤ x ≤ maxs[i]
          push!(yd,spl[i](x))
        else
          push!(yd,NaN)
        end
      end
    end
    # Add a column with completed y data of current DataFrame to the output matrix
    ydata = hcat(ydata,yd)
  end

  # Return unified y data
  return ydata
end #function get_ydata


"""
get_ScalingFactor(lower_xdata, common_xdata, upper_xdata,
                  lower_ydata, common_ydata, upper_ydata, mi, Mi)

From the unified x and y data below, inside, and above the common range and
the indices of ending vectors outside the common range in the unified x data
(`mi` and `Mi`), calculate and return scaling factors to allow a smooth contineous
mean at the premature ends of vectors.
"""
function get_ScalingFactor(lower_xdata, common_xdata, upper_xdata,
                           lower_ydata, common_ydata, upper_ydata, mi, Mi)

  ### Get scaling factors for data below the common x range
  flow = ones(lower_xdata); fhigh = ones(upper_xdata)
  if isempty(mi) && length(lower_ydata) > 0
    # Calculate scaling factors, if no premature ending vectors are found
    # outside the common range
    flow[1:end] = mean(common_ydata[1,:])/
                  mean(common_ydata[1,:][!isnan.(lower_ydata[end,:])])
  elseif length(lower_ydata) > 0
    # Calculate scaling factors, if premature ending vectors outside the
    # common range are found
    y = lower_ydata[mi[1],:]
    # Scaling factors from common range to first ending vectors
    flow[1:mi[1]-1] = mean(y[!isnan.(y)])/mean(y[!isnan.(lower_ydata[mi[1]-1,:])])
    # Scaling vectors for further ending vectors
    for i = 2:length(mi)
      y = lower_ydata[mi[i],:]
      flow[mi[i-1]:mi[i]-1] = mean(y[!isnan.(y)])/
                              mean(y[!isnan.(lower_ydata[mi[i]-1,:])])
    end
    # Scaling factors for last prematurely ending vector to end of data range
    flow[mi[end]:end] = mean(common_ydata[1,:])/
                        mean(common_ydata[1,:][!isnan.(lower_ydata[end,:])])
  end

  ### Get scaling factors for data above the common x range
  if isempty(Mi) && length(upper_ydata) > 0
    # Calculate scaling factors, if no premature ending vectors are found
    # outside the common range
    fhigh[1:end] = mean(common_ydata[end,:])/
                   mean(common_ydata[end,:][!isnan.(upper_ydata[1,:])])
  elseif length(upper_ydata) > 0
    # Calculate scaling factors, if premature ending vectors outside the
    # common range are found
    # Scaling factors from common range to first ending vectors
    fhigh[1:Mi[1]-1] = mean(common_ydata[end,:])/
                       mean(common_ydata[end,:][!isnan.(upper_ydata[1,:])])
    # Scaling vectors for further ending vectors
    for i = 2:length(Mi)
      y = upper_ydata[Mi[i],:]
      fhigh[Mi[i-1]:Mi[i]-1] = mean(y[!isnan.(y)])/
                               mean(y[!isnan.(upper_ydata[Mi[i]+1,:])])
    end
    # Scaling vectors for further ending vectors
    y = upper_ydata[Mi[end],:]
    fhigh[Mi[end]:end] = mean(y[!isnan.(y)])/
    mean(y[!isnan.(upper_ydata[Mi[end]+1,:])])
  end

  # Return scaling factors for mean below and above the common range
  return flow, fhigh
end #function get_ScalingFactor


"""
    find_boundaries(mins,maxs,lower_xdata,upper_xdata)

From the minima (`mins`) and maxima (`maxs`) in the x data of each DataFrame,
and the `lower_xdata` and `upper_xdata` below and above the common range, return
all indices in `lower_xdata` and `upper_xdata`, where vectors are prematurely
ending outside the common range.
"""
function find_boundaries(mins,maxs,lower_xdata,upper_xdata)
  # Initialise
  mi = Int64[]; Mi = Int64[]
  # Loop over extrema
  for m = 1:length(mins)
    # Save vectors starting later the absolute minimum,
    # but before common range
    if all(float(mins[m]).!=extrema(mins))
      push!(mi,find(lower_xdata.==mins[m])[1])
    end
    # Save prematurely ending vectors outside common range
    if all(float(maxs[m]).!=extrema(maxs))
      push!(Mi,find(upper_xdata.==maxs[m])[1])
    end
  end
  # Unify data for more than one ending vector at the same point
  mi = sort(unique(mi)); Mi = sort(unique(Mi))

  # Return indices
  return mi, Mi
end


"""
    calc_MeanSum(Ydata,SF)

From the unified `Ydata` of each DataFrame and the corresponding scaling factors
`SF` for the mean outside the common range, calculate the `mean` and `sum` of
all `Ydata`.
"""
function calc_MeanSum(Ydata,SF)
  Ymean = Float64[]; Ysum = Float64[]
  for i = 1:length(SF)
    push!(Ymean,mean(Ydata[i,:][!isnan.(Ydata[i,:])])⋅SF[i])
    push!(Ysum,sum(Ydata[i,:][!isnan.(Ydata[i,:])]))
  end

  return Ymean, Ysum
end #function calc_MeanSum

end #module VectorMaths
