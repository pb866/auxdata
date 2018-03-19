module numerics

export int, dig2int


"""
    int(n)

Converts a `String` or a `Float` to `Int64`. The `Float` must not contain after
period digits and the `String` must contain a valid integer.
"""
function int(n)
  return convert(Int64,float(n))
end #function int


"""
    dig2int(dig::Array{Int64,1})

Transform an array of integer digits derived from digits(itg) back to an integer.
"""
function dig2int(dig::Array{Int64,1},base::Int64=10)

  # Initialise
  itg = 0
  # Loop over digits (starting with first digit at the end of the array)
  for i = length(dig):-1:1
    # Add digit to overall integer as multiples of 10^(n-1), n = position in integer
    itg += dig[i]*base^(i-1)
  end
  return itg
end #function dig2int

end #module numerics
