__precompile__()


"""
# Data Analysis and Manipulation tools
"""
module VectorMaths

# Load packages
try using DataFrames
catch
  Pkg.add("DataFrames")
  using DataFrames
end
try Dierckx
catch
  Pkg.add("Dierckx")
  using Dierckx
end

export DFvector_mean


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

end #module VectorMaths
