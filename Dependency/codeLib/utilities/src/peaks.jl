export peaks
# 1. The Scalar Function (Base definition)
function peaks_loc(x::Real, y::Real)
    z =  3 * (1 - x)^2 * exp(-(x^2) - (y+1)^2) - 
         10 * (x/5 - x^3 - y^5) * exp(-x^2 - y^2) - 
         1/3 * exp(-(x+1)^2 - y^2)
    return z
end

# 2. Vector Input Method (Handles custom x and y vectors)
#    This effectively does the 'meshgrid' calculation via broadcasting.
function peaks(x::AbstractVector, y::AbstractVector)
    # Broadcast scalar function over x' (1xN) and y (Mx1)
    # Result is an MxN matrix
    return peaks_loc.(x, y)
end

# 3. Convenience Method for Dimensions (nx, ny)
#    Generates the standard -3 to 3 range with defined lengths.
function peaks(nx::Int, ny::Int)
    x = range(-3, 3, length=nx)
    y = range(-3, 3, length=ny)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return peaks_loc.(X, Y)
end

# 4. Fallback for Square Grid (n)
peaks(n::Int=49) = peaks(n, n)