using DifferentialEquations
using Plots
using CSV
using FFTW
using StatsBase
using DataFrames
using Base.Threads

filepostion = "/home/gnugnu/문서/Pendulum-experiment-analysis"
print("filename: ")
filename = readline()
data = CSV.read(filepostion*"/data/"*filename*".csv") 
light = data[!,:2]  
n=length(light)
println("n: "*string(n))

newlight = Float64[]
for i in 1:n
    push!(newlight, light[i] - minimum(light))
end

l = 0.2                             # length [m]
R = 0.012/2                         # radius of ball 
m = 7.1206e-3                       # mass of ball 
r = 0.007/2                         # radius of light sensor
g = 9.81                            # gravitational acceleration [m/s²]
dt = 2e-3
print("θ₀[rad]: ")
theta0 = parse(Float64, readline())
θ₀ = atan(theta0)                    # initial angular deflection [rad]
ω₀ = 0.0                           # initial angular velocity [rad/s]
u₀ = [θ₀, ω₀]                       # initial state vector
tspan = (0.0, data[n,:1])        # time interval
tstop = 0.0:dt:data[n,:1]

function f(theta) # i have to fix because ball and light's radius is different 
    a = l*sin(abs(theta/2))
    if(a<r) s = a*sqrt(r^2-a^2) + r^2*atan(a/(r^2-a^2))
        else s = pi*r^2
    end
    return s 
end

function f2(theta)
  s = maximum(newlight);
  a = l*sin(theta)
  b = l - l*cos(theta)
  x1 = sqrt(R^2-r^2) - a
  x2 = -sqrt(R^2-r^2) - a

  ans = 0;
  if(x1<=-r) 
    return s
  end
  if(-r<=x1 && x1<=r) 
    ans += (r-x1)/(2r)
  end
  if(-r<=x2 && x2<=r) 
    ans += (r+x2)/(2r)
  end
  if(x2>=r) 
    return s;
  end
  return ans*s
end

function pendulum(du,u,p,t)
    du[1] = u[2]                    # θ'(t) = ω(t)
    du[2] = -g/l*sin(u[1]) - p(t)*(u[2])/m  # ω'(t) = -g/l sin θ(t) - b ω(t) / m
end

Fexp = fft(newlight)
absFexp = Array{Float64}(undef, length(Fexp))
for i in 1:length(Fexp)
    absFexp[i] = abs(Fexp[i])
end

println("start fftbm")
num = 0.0:0.5:20
AA = Array{Float64}(undef, length(num))
bb = Array{Float64}(undef, length(num))
for i in 1:length(num)
    b = t -> 6*pi*1e-3*r+num[i]*1e-5
    bb[i] = b(0)
    prob = ODEProblem(pendulum,u₀,tspan,b)
    alg = RK4()
    sol = solve(prob,alg,dt = 2e-3,adaptive=false)
    result = Array{Float64}(undef, length(sol.t))
    for j in 1:length(sol.t)
        result[j] = f2(sol[1,j])
    end
    Fcal = fft(result[1:n])
    absdeltaF = Array{Float64}(undef, length(Fexp))
    for j in 1:length(Fexp)
        absdeltaF = abs(Fexp[j]-Fcal[j])
    end
    x = sum(absdeltaF)/sum(absFexp)
    global AA[i] = x
end
plot(bb, AA, title = "FFTBM", label = "AA")
xlabel!("b (1/s)")
ylabel!("AA")
name = filepostion*"/result/"*filename*"/"*filename*"_AA.png"
savefig(name)
println("end fftbm")

b = t -> 6*pi*1e-3*r + num[findfirst(x->x==minimum(AA),AA)]*1e-5
bstring = string(b(0))
prob = ODEProblem(pendulum,u₀,tspan,b)
alg = RK4()
sol = solve(prob,alg,tstops = tstop ,adaptive=false)

result = Float64[]
for i in 1:length(newlight)
    push!(result, f2(sol[1,i]))
end

c = crosscor(newlight, result)
cmax = maximum(c)
c = c/cmax
t = findfirst(x->x==1.0, c)
println("t: "*string(t))

function graph(start, d)

    plot(sol.t[start:start+d], newlight[start:start+d], label = "data", title = "simulation at b = "*bstring, linealpha = 1.0)
    plot!(sol.t[start:start+d],result[start-t:start+d-t],label = "simulation", linealpha = 0.7)
    xlabel!("Time(s)")
    ylabel!("Light Intensity(%)")
    savefig(filepostion*"/result/"*filename*"/"*filename*" SIM("*string(start*dt)*", "*string((start+d)*dt)*").png")
end

print("how many?")
number = parse(Int64, readline())
for i in 1:number
print("start")
start = parse(Int64, readline())
print("d")
d = parse(Int64, readline())
graph(start, d)
end

df = DataFrame(Time = sol.t[t+1:n], Data = newlight[t+1:n], Simulation = result[1:n-t])
print(df)
CSV.write(filepostion*"/result/"*filename*"/"*filename*"_result.csv", df)

