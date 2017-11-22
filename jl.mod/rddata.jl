__precompile__()


"""
# Module _rddata_

Read in data from TUV (version 5.2 format) output file and save
χ-dependent _j_ values to dataframe.

# Functions

- read_j (public)
- read_head (private)
- read_data (private)
- read_between (private)
"""
module rddata
using DataFrames
export read_j


"""
    read_j(<TUV 5.2 input file>)

Read in data from TUV (version 5.2 format) output file and save
χ-dependent _j_ values to dataframe.
"""
function read_j(ifile)

  # Read reactions and j values from input file
  open(ifile,"r") do f
    lines = readlines(f)
    jvals = read_head(lines)
    jvals = read_data(lines,jvals)
    return jvals
  end
end #function read_j


"""
    read_head(<TUV 5.2 input file>)

Initialise dataframe with heaer.
- First column: solar zenith angles (sza, χ)
- Further columns: χ-dependent _j_ values
"""
function read_head(lines)

  # Retrieve lines from TUV file between key words omitting labels before equal sign
  rxns = read_between(lines, "Photolysis rate coefficients, s-1", "values at z", " = ", 2)

  # Initialise dataframe with sza and reactions columns
  jvals = DataFrame()
  header = Symbol["sza"]
  for j in rxns  push!(header,j)  end
  for name in header  jvals[name] = Float64[]  end

  # Return dataframe
  return jvals
end # function read_head


"""
    read_data(lines, df)

From the lines (lines) in the TUV file retrieve χ-dependent _j_ values
and append dataframe (df) initialised by function _read_head_.
"""
function read_data(lines,jvals)

  # Array with data lines as strings
  rawdata = read_between(lines,"sza, deg.", "---")
  # Initialise j data as floats
  jdata = []
  # Loop over raw data and split into matrix of floats
  for i in 1:length(rawdata)
    data = split(rawdata[i])
    dataline = map(x->parse(Float64,x),data)
    push!(jdata,dataline)
  end
  # Append dataframe of j values
  for j in jdata  push!(jvals,j)  end

  # Return completed dataframe
  return jvals

end #function read_data


"""
    read_between(lines, lstart, lend[[, Sep], col])

Return the lines of strings from all lines of a file (lines) in-between the lines
containing the key phrases 'lstart' and 'lend'.

If an optional separator (Sep) and a column number (col) is defined, lines will be
split by Sep and only column col is returned.
"""
function read_between(lines, lstart::String, lend::String, Sep::String="", col::Int64=0)

  # Find first line with with reaction definitions in TUV input file
  start = Int64 # Initialise index of starting line to use outside of loop
  for i = 1:length(lines)
    if contains(lines[i], lstart)
      start = i + 1
      break
    end
  end

  # Create array with reactions from definition list in TUV input file
  selected_lines = Array{String}(0)
  for i = start:length(lines)
    if contains(lines[i], lend)  break  end
    if col ≠ 0  lines[i] = split(lines[i],Sep)[col]  end
    push!(selected_lines,strip(lines[i]))
  end
  return selected_lines
end #function read_between

end # module rddata
