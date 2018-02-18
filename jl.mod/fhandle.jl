__precompile__()


"""
# Module _fhandle_

Routines for file handling, reading/storing input data and writing output data
to text files.

# Functions

- get_fnames
- wrt_params
- test_dir
- test_file
- rdfil
- rdinp
"""
module fhandle
try using DataFrames
catch
  Pkg.add("DataFrames")
  using DataFrames
end
try using Juno: input
catch
  Pkg.add("Juno")
  using Juno: input
end

export get_fnames,
       wrt_params,
       test_dir,
       test_file,
       rdfil,
       rdinp


"""
    get_fnames(scen)

From the scenario name (scen) used in TUV derive names for the TUV output file
and the output folder.
"""
function get_fnames(scen)

  # Test script argument and ask for input file, if missing
  if scen == ""
    println("Enter file with TUV scenario name")
    scen = input("(Name of output file without '.txt'): ")
  else
    scen = ARGS[1]
  end

  # Create output folder
  iofolder = "./params_"*strip(scen)
  try mkdir(iofolder)
  catch
    print("\033[95mFolder '$iofolder' already exists! Overwrite ")
    confirm = input("(\033[4mY\033[0m\033[95mes/\033[4mN\033[0m\033[95mo)?\033[0m ")
    if lowercase(confirm) == "yes" || lowercase(confirm) == "y"
      cd(iofolder); files = readdir(); for file in files  rm(file)  end; cd("..")
    else println("Programme aborted. Exit code '98'."); exit(98)
    end
  end

  # Define TUV file
  ifile = scen*".txt"
  if isfile(ifile) == false
    println("$ifile does not exist. Script terminated. Exit code: 99."); exit(99)
  end

  # return file and folder names
  return ifile, iofolder
end #function get_fnames


"""
    wrt_params(rxn, fit, sigma, rmse, R2, iofolder, time)

For each reaction (rxn), from param in fit and sigma, print parameters
and 95% confidence together with RMSE and R^2 to file 'parameters.dat'
in the designated output folder (iofolder).

Print the time of creation (time) to the output file.
"""
function wrt_params(rxn, fit, sigma, rmse, R2, iofolder, time)

  #transform dataframe column symbols to strings
  rxn = convert.(String,rxn)

  # Open output file
  open("$iofolder/parameters.dat","w") do f
    # Print header
    println(f,
    "Parameters and statistical data for parameterisation of photolysis processes")
    println(f,
    "in the Master Chemical Mechanism (MCM; http://mcm.leeds.ac.uk/MCMv3.3.1/) using:")
    println(f, "\nj / s-1 = l·(cos(x))^m·exp(-n·sec(x))")
    println(f, "\n\ncreated $(Dates.format(time,"dd.mm.yyyy, HH:MM:SS"))")
    println(f,"\n                P a r a m e t e r s               S t a t i s t i c s")
    println(f,"    l / s-1              m              n         RMSE / s-1    R^2      Reaction")

    # Loop over reactions
    for i = 1:length(fit)
      # Determine order of magnitude of parameter l
      p = floor(log10(fit[i].param[1]))
      # Print parameters, statistical data, and reaction label to output file
      @printf(f,"(%.3f±%.3f)E%d    %.3f±%.3f    %.3f±%.3f    %.3e    %.4f    %s\n",
      fit[i].param[1]/10^p,sigma[i][1]/10^p,p,
      fit[i].param[2],sigma[i][2],fit[i].param[3],sigma[i][3],rmse[i],R2[i],rxn[i+1])
    end
  end
end # function wrt_params


"""
    test_dir(ifile::String;default_dir::String="./")

Check for existance of the directory specied in `ifile` after the removal of the
file name. If directory doesn't exist, ask to create directory or specify a new
path. If `default_dir` is specified, rdinp will look for `ifile` in this
directory, if `ifile` does not include a folder path.
"""
function test_dir(ifile::String; default_dir::String="./")

  # Add default directory, if folder path in file name is missing
  dir = dirname(ifile)
  if dir == ""
    dir = default_dir
    ifile = normpath(joinpath(dir,basename(ifile)))
  end

  # Test directory path and warn, if non-existant, either create path (opt. 1)
  # redefine output file + path (opt. 2) or exit script (general case)
  if !isdir(dir)
    println("Directory for specified file does not exist!")
    println("Create directory [1] or specify new (folder +) file [2] or quit [<ENTER>]?")
    decision = input()
    if decision == "1"
      mkpath(dir)
    elseif decision == "2"
      ifile = input("Enter file name: ")
      ifile = test_dir(ifile,default_dir=default_dir)
    else
      exit()
    end
  end

  # Return (revised) file name (+ path)
  return ifile
end #function test_file


"""
    test_file(ifile::String;default_dir::String="./")

Check for existance of ifile. If file doesn't exist, ask for a file name until
file is found. If `default_dir` is specified, rdinp will look for `ifile` in this
directory, if `ifile` does not include a folder path.
"""
function test_file(ifile;default_dir::String="./")

  # Add default directory, if folder path in file name is missing
  fname = basename(ifile); dir = dirname(ifile)
  if dir == ""  dir = default_dir  end
  ifile = joinpath(dir,fname)

  # Test existance of file or ask for user input until file is found
  while !isfile(ifile)
    println("File $ifile does not exist!")
    ifile = input("Enter file (or press <ENTER> to quit): ")
    if ifile==""  exit()  end
  end

  return ifile
end #function test_file


"""
    rdfil(ifile::String,rmhead::Int64=0)

Read content from ifile return an array with its lines.

If specified, the first `rmhead` lines are omitted.
"""
function rdfil(ifile::String,rmhead::Int64=0)

  lines =[]
  # Read lines from file
  open(ifile,"r") do f
    lines = readlines(f)
    # Delete header
    for i = 1:rmhead
      shift!(lines)
    end
  end

  # Return array with lines
  return lines
end #function rdinp


"""
    rdinp(ifile::String,rmhead::Int64=0;default_dir::String="./")

Check for existance of file and return array with its lines.

The first `rmhead` lines are omitted (none line by default).
If `default_dir` is specified, rdinp will look for `ifile` in this directory,
if `ifile` does not include a folder path.
"""
function rdinp(ifile::String,rmhead::Int64=0;default_dir::String="./")

  # Add default directory, if folder path in file name is missing
  fname = basename(ifile); dir = dirname(ifile)
  if dir == ""  dir = default_dir  end
  ifile = normpath(joinpath(dir,fname))

  # Check existance of database file
  while !isfile(ifile)
    println("File $ifile does not exist!")
    ifile = input("Enter file (or press <ENTER> to quit): ")
    if ifile==""  exit()  end
  end

  # Read lines from database file and delete header
  lines =[]
  open(ifile,"r+") do f
    lines = readlines(f)
    for i = 1:rmhead
      shift!(lines)
    end
  end

  # Return array with lines
  return lines
end #function rdinp
end # module fhandle
