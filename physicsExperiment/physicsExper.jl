using DifferentialEquations
using Plots
using CSV

data = CSV.read("/home/gnugnu/문서/physicsExperiment/num4.csv")
light = data[!,:2]  
n=length(light)
min = max = light[1]
for i in 1:n
    global  min, max    
    if(min > light[i]) min = light[i]
    end
    if(max < light[i]) max = light[i]
    end
end

newlight = Float64[]
for i in 1:n
    push!(newlight, light[i] - min)
end

l = 0.2                             # length [m]
R = 0.012/2                             # radius of ball 
m = 7.1206e-3
r = 0.007/2                             #radius of detector
g = 9.81                            # gravitational acceleration [m/s²]
b = 6*pi*1e-3*r+2.5e-5


function pendulum!(du,u,p,t)
    du[1] = u[2]                    # θ'(t) = ω(t)
    du[2] = -g/l*sin(u[1]) - b*(u[2])/m  # ω'(t) = -g/l sin θ(t) - b ω(t) / m
end

θ₀ = 0.049958395721942765        # initial angular deflection [rad]
ω₀ = 0.0                           # initial angular velocity [rad/s]
u₀ = [θ₀, ω₀]                       # initial state vector
tspan = (0.0, data[n,:1])                 # time interval


prob = ODEProblem(pendulum!,u₀,tspan)
alg = RK4()
sol = solve(prob,alg,dt = 2e-3,adaptive=false)

function f(theta) # i have to fix because ball and detector's radius is different 
    a = l*sin(abs(theta/2))
    if(a<r) s = a*sqrt(r^2-a^2) + r^2*atan(a/(r^2-a^2))
        else s = pi*r^2
    end
    return s 
end

function f2(theta)
  s = max-min;
  a = l*sin(theta)
  b = l - l*cos(theta)
  x1 = sqrt(R^2-r^2) - a
  x2 = -sqrt(R^2-r^2) - a

  ans = 0;
  if(x1<=-r) return s
  end
  if(-r<=x1 && x1<=r) ans += (r-x1)/(2r)
  end
  if(-r<=x2 && x2<=r) ans += (r+x2)/(2r)
  end
  if(x2>=r) return s;
  end
  return ans*s
end


result = Float64[]
for i in 1:length(sol.t)
    push!(result, f2(sol[1,i]))
end

plot(sol.t,result)
plot!(sol.t, newlight)
savefig("/home/gnugnu/문서/physicsExperiment/wa2.png")

