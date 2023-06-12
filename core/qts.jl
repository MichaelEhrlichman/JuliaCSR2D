using IntelVectorMath

function samplepoint_will(t::Float64)
    #println("t: ", t)
    sinht = sinh(t)
    ϕ = tanh(sinht*π/2)
    ϕ′ = (cosh(t)*π/2)/cosh(sinht*π/2)^2
    return ϕ, ϕ′
end

function samplepoint_fast(t::Float64)
  et = exp(t)
  et2 = et^2
  pt_sinht = pi/4*(et2-1)/et  # pi/2*sinh(t)
  pt_cosht = pi/4*(et2+1)/et  # pi/2*cosh(t)

  bet=exp(pt_sinht)
  bet2=bet^2
  phi = (bet2-1)/(bet2+1)
  phip = pt_cosht/( (bet2+1)/2/bet )^2
  return phi, phip
end


function samplepoint_V(t::Vector{Float64})
  et = IVM.exp(t)
  et2 = @. et^2
  pt_sinht = @. pi/4*(et2-1.0)/et
  pt_cosht = @. pi/4*(et2+1.0)/et

  bet = IVM.exp(pt_sinht)
  bet2 = @. bet^2
  phi = @. (bet2-1.0)/(bet2+1.0)
  phip = @. pt_cosht/( (bet2+1.0)/2.0/bet )^2
  return phi, phip
end

using LinearAlgebra: norm
function estimate_error_will(prevI, I)
    ε = eps(Float64)
    M = 20
    return M*norm(I - prevI)^2 + norm(I)*ε
end

# OBSOLETE
function qts_will_2(f::Function, h0::Float64, maxlevel::Integer=6)

    #@assert maxlevel > 0
    #@assert h0 > 0
    
    T = Float64
    
    x0, w0 = samplepoint_will(0.0)  #origin
    Σ = f(x0)*w0
   
    k = 1
    while true
        #t = k*h0/2^maxlevel
        xk, wk = samplepoint_will(k*h0/2^maxlevel)  
        1 - xk ≤ eps(T) && break   # xk is too close to 1, the upper bound of integral
        wk ≤ floatmin(T) && break  # wk is too small, series trucated
        
        Σ += (f(xk) + f(-xk)) * wk
        k += 1    # step is either 1 (for level = 0) or 2

    end    
    h = h0/2^maxlevel
    I = h*Σ
    
    #h = h0/2^maxlevel
    #return h*Σ
    return I
    #return I, E

end

function qts_V(f::Function, h0::Float64, maxlevel::Integer=6)

    T = Float64
    ut_xk = 3.1540  # asinh(2/pi*atanh(1-eps(Float64)))
    ut_wk = 6.1227  # found numerically with mathematica for φ'(t)≈floatmin(Float64)
    
    x0,w0 = samplepoint_V([0.0])
    Σ = f(x0[1])*w0[1]
    
    I=0
    k=1
    for level in 0:maxlevel
        h = h0/2^level
        Nt = floor(min(ut_xk, ut_wk) / h)
        tvec = [ h+h*n for n=0:k:Nt]
        k=2 # step is either 1 (for level = 0) or 2
        x, w = samplepoint_V(tvec)
        for (xk,wk) in zip(x,w)
            Σ += (f(xk) + f(-xk)) * wk
        end
        if level == 0
            I=h*Σ
        else
            prevI = I
            I=h*Σ
            E=estimate_error_will(prevI,I)
            tol = norm(I)*sqrt(eps(T))
            !(E>tol) && break
        end
    end
    return I
end

function qts_will(f::Function, h0::Float64, maxlevel::Integer=6)

    #@assert maxlevel > 0
    #@assert h0 > 0
    T = Float64
    x0, w0 = samplepoint_will(0.0)  #origin
    Σ = f(x0)*w0
   
    k = 1
    while true
        t = k*h0
        xk, wk = samplepoint_will(t)
        1 - xk ≤ eps(T) && break   # xk is too close to 1, the upper bound of integral
        wk ≤ floatmin(T) && break  # wk is too small, series trucated
        
        Σ += (f(xk) + f(-xk)) * wk
        k += 1    # step is either 1 (for level = 0) or 2

    end    
    I = h0*Σ
    E = zero(eltype(I))
    for level in 1:maxlevel
        k = 1
        h = h0/2^level
        while true
            t = k*h
            xk, wk = samplepoint_will(t)
            
            1 - xk ≤ eps(T) && break   # xk is too close to 1, the upper bound of integral
            wk ≤ floatmin(T) && break  # wk is too small, series trucated
            Σ += (f(xk) + f(-xk)) * wk
            k += 2     # step is either 1 (for level = 0) or 2

        end     
        prevI = I
        I = h*Σ
        E = estimate_error_will(prevI, I)
        ###tol = max(norm(I)*rtol, atol)
        tol = norm(I)*sqrt(eps(T))
        !(E > tol) && break
    end
    #h = h0/2^maxlevel
    #return h*Σ
    return I
    #return I, E
end


# In contrast to qts_will which integrates from -1 to 1
# QTS_will integrates from a to b
function QTS_will(f::Function, a::Float64, b::Float64, M::Int)
    s = (b + a)/2
    t = (b - a)/2

    #print("t:", t)
    #I, E = q.qts(u -> f(s + t*u); atol=atol/t, rtol=rtol)
    #I*t, E*t
    
    #h0 = 1.0
    #maxlevel = 12 # A large maxlevel may be slow if the integral does not converge...
    I = qts_will(u-> f(s+t*u), 1.0, M)
    ##I = qts_V(u-> f(s+t*u), 1.0, M)
    
    return I*t   

end 