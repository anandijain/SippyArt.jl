function lotka_volterra(du,u,p,t)
  du[1] = p[1]*u[1] - p[2]*u[1]*u[2] 
  du[2] = -p[3]*u[2] + p[4]*u[1]*u[2] 
end

function lorenz(du, u, p, t)
    du[1] = p[1] * (u[2] - u[1])
    du[2] = u[1] * (p[2] - u[3]) - u[2]
    du[3] = u[1] * u[2] - p[3] * u[3]
end

function duffing(du,u,p,t)
 du[1] = u[2]
 du[2] = u[1] - u[1]^3 -p[1]*u[2] + p[2]*cos(p[3]*t)
end

function cyclically_symmetric(du,u,p,t)
 du[1] = sin(u[2]) - p[1]*u[1]
 du[2] = sin(u[3]) - p[1]*u[2]
 du[3] = sin(u[1]) - p[1]*u[1]
end

function rossler(du,u,p,t)
 du[1] = -u[2] - u[3]
 du[2] = u[1] + p[1]*u[2]
 du[3] = p[2] * u[3]*(u[1] - p[3])
end

function hindmarsh_rose(du,u,p,t)
 du[1] = u[2] - p[1]*u[1]^3+p[2]*u[1]^2-u[3]+p[3]
 du[2] = p[4]-p[5]*u[1]^2 - u[2]
 du[3] = p[6]*(p[7]*(u[1] - p[8]) - u[3])
end

function aizawa(du, u, p, t)
    du[1] = (u[3]-p[2]) * u[1] - p[4]*u[2]
    du[2] = p[4] * u[1] + (u[3]-p[2]) * u[2]
    du[3] = p[3] + p[1]*u[3] - u[3]^3 / 3 - (u[1]^2 + u[2]^2)*(1 + p[5]*u[3]) + p[6]*u[3]*u[1]^3
end


function chen(du, u, p, t)
  du[1] = p[1]*u[1]-u[1]*u[2]
  du[2] = p[2]*u[2]+u[1]*u[3]
  du[3] = p[3]*u[3]+u[1]*u[2]/3
end

# function halvorsen(du, u, p, t)
#   du[1] = -p[1]*u[1]
#   du[2] =
#   du[3] =
# end

# function sprott(du, u, p, t)
#   du[1] = 
#   du[2] =
#   du[3] =
# end




