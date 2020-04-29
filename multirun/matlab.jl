using MATLAB

# example array
array = [1,2,3,4,5];

# to call a Julia variable or array directly from MATLAB (the conversion happens on the fly), just place $ in front

# MATLAB script section, use " " for one line, """ """ for multiline
mat"""
double_array = 2 * $array;
plot(double_array)
"""

# to retrieve data from MATLAB, use @mget
jarray = @mget double_array;

println(jarray)

# Only use MATLAB integration for math and plots, printing freaks it out
# Ignore the documentation for matlab.jl also