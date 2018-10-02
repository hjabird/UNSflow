#Function for estimating a problem's time step
function update_a2a3adot(surf::TwoDSurf,dt)
    for ia = 2:3
        surf.aterm[ia] = simpleTrapz(surf.downwash.*cos.(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_a2a3adot(surf::TwoDSurf2DOF,dt)
    for ia = 2:3
        surf.aterm[ia] = simpleTrapz(surf.downwash.*cos.(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_adot(surf::TwoDSurf,dt)
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

function update_adot(surf::TwoDSurf2DOF,dt)
    surf.a0dot[1] = (surf.a0[1] - surf.a0prev[1])/dt
    for ia = 1:3
        surf.adot[ia] = (surf.aterm[ia]-surf.aprev[ia])/dt
    end
    return surf
end

# Function for updating the induced velocities
function update_indbound(surf::TwoDSurf, curfield::TwoDFlowField)
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev; curfield.extv], surf.bnd_x, surf.bnd_z)
    return surf
end

function update_indbound(surf::TwoDSurf2DOF, curfield::TwoDFlowField)
    surf.uind[1:surf.ndiv], surf.wind[1:surf.ndiv] = ind_vel([curfield.tev; curfield.lev; curfield.extv], surf.bnd_x, surf.bnd_z)
    return surf
end

# Function for updating the downwash
function update_downwash(surf::TwoDSurf, vels::Vector{Float64})
    for ib = 1:surf.ndiv
        surf.downwash[ib] = -(surf.kinem.u + vels[1])*sin(surf.kinem.alpha) - surf.uind[ib]*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) - surf.wind[ib]*cos(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[ib] - surf.pvt*surf.c) + surf.cam_slope[ib]*(surf.uind[ib]*cos(surf.kinem.alpha) + (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
    end
    return surf
end

function update_downwash(surf::TwoDSurf2DOF, vels::Vector{Float64})
    for ib = 1:surf.ndiv
        surf.downwash[ib] = -(surf.kinem.u + vels[1])*sin(surf.kinem.alpha) - surf.uind[ib]*sin(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*cos(surf.kinem.alpha) - surf.wind[ib]*cos(surf.kinem.alpha) - surf.kinem.alphadot*(surf.x[ib] - surf.pvt*surf.c) + surf.cam_slope[ib]*(surf.uind[ib]*cos(surf.kinem.alpha) + (surf.kinem.u + vels[1])*cos(surf.kinem.alpha) + (surf.kinem.hdot - vels[2])*sin(surf.kinem.alpha) - surf.wind[ib]*sin(surf.kinem.alpha))
    end
    return surf
end

# Function for a_0 and a_1 fourier coefficients
function update_a0anda1(surf::TwoDSurf)
    surf.a0[1] = simpleTrapz(surf.downwash,surf.theta)
    surf.aterm[1] = simpleTrapz(surf.downwash.*cos.(surf.theta),surf.theta)
    surf.a0[1] = -surf.a0[1]/(surf.uref*pi)
    surf.aterm[1] = 2.*surf.aterm[1]/(surf.uref*pi)
    return surf
end

function update_a0anda1(surf::TwoDSurf2DOF)
    surf.a0[1] = simpleTrapz(surf.downwash,surf.theta)
    surf.aterm[1] = simpleTrapz(surf.downwash.*cos.(surf.theta),surf.theta)
    surf.a0[1] = -surf.a0[1]/(surf.uref*pi)
    surf.aterm[1] = 2.*surf.aterm[1]/(surf.uref*pi)
    return surf
end

# Function for calculating the fourier coefficients a_2 upwards to a_n
function update_a2toan(surf::TwoDSurf)
    for ia = 2:surf.naterm
        surf.aterm[ia] = simpleTrapz(surf.downwash.*cos.(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    return surf
end

function update_a2toan(surf::TwoDSurf2DOF)
    for ia = 2:surf.naterm
        surf.aterm[ia] = simpleTrapz(surf.downwash.*cos.(ia*surf.theta),surf.theta)
        surf.aterm[ia] = 2.*surf.aterm[ia]/(surf.uref*pi)
    end
    return surf
end


#Function to update the external flowfield
function update_externalvel(curfield::TwoDFlowField, t)
    if (typeof(curfield.velX) == CosDef)
        curfield.u[1] = curfield.velX(t)
        curfield.w[1] = curfield.velZ(t)
    elseif (typeof(curfield.velX) == SinDef)
        curfield.u[1] = curfield.velX(t)
        curfield.w[1] = curfield.velZ(t)
    elseif (typeof(curfield.velX) == ConstDef)
        curfield.u[1] = curfield.velX(t)
        curfield.w[1] = curfield.velZ(t)
    end
    if (typeof(curfield.velX) == StepGustDef)
        curfield.u[1] = curfield.velX(t)
    end
    if (typeof(curfield.velZ) == StepGustDef)
        curfield.w[1] = curfield.velZ(t)
    end
end



# Function updating the dimensional kinematic parameters
function update_kinem(surf::TwoDSurf, t)

    # Pitch kinematics
    if (typeof(surf.kindef.alpha) == EldUpDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == EldUptstartDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == EldRampReturnDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == ConstDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = 0.
    elseif (typeof(surf.kindef.alpha) == SinDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == CosDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    elseif (typeof(surf.kindef.alpha) == VAWTalphaDef)
        surf.kinem.alpha = surf.kindef.alpha(t)
        surf.kinem.alphadot = ForwardDiff.derivative(surf.kindef.alpha,t)*surf.uref/surf.c
    end
    # ---------------------------------------------------------------------------------------------

    # Plunge kinematics
    if (typeof(surf.kindef.h) == EldUpDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUptstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpIntDef)
            surf.kinem.h = surf.kindef.h(t)*surf.c
            surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldUpInttstartDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == EldRampReturnDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == ConstDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = 0.
    elseif (typeof(surf.kindef.h) == SinDef)
      surf.kinem.h = surf.kindef.h(t)*surf.c
      surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == CosDef)
      surf.kinem.h = surf.kindef.h(t)*surf.c
      surf.kinem.hdot =  ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
    elseif (typeof(surf.kindef.h) == VAWThDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
        end
    # ---------------------------------------------------------------------------------------------

    # Forward velocity
    if (typeof(surf.kindef.u) == EldUpDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == EldRampReturnDef)
        surf.kinem.u, surf.kinem.udot = surf.kindef.u(t)
        surf.kinem.u = surf.kinem.u*surf.uref
        surf.kinem.udot = surf.kinem.udot*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == ConstDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = 0.
    elseif (typeof(surf.kindef.h) == VAWThDef)
        surf.kinem.h = surf.kindef.h(t)*surf.c
        surf.kinem.hdot = ForwardDiff.derivative(surf.kindef.h,t)*surf.uref
        end
    # ---------------------------------------------------------------------------------------------

    # Forward velocity
    if (typeof(surf.kindef.u) == EldUpDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == EldRampReturnDef)
        surf.kinem.u, surf.kinem.udot = surf.kindef.u(t)
        surf.kinem.u = surf.kinem.u*surf.uref
        surf.kinem.udot = surf.kinem.udot*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == ConstDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = 0.
    elseif (typeof(surf.kindef.u) == SinDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == CosDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == LinearDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    elseif (typeof(surf.kindef.u) == VAWTuDef)
        surf.kinem.u = surf.kindef.u(t)*surf.uref
        surf.kinem.udot = ForwardDiff.derivative(surf.kindef.u,t)*surf.uref*surf.uref/surf.c
    end
    # ---------------------------------------------------------------------------------------------
    return surf
end

# ---------------------------------------------------------------------------------------------
# Updates the bound vorticity distribution: eqn (2.1) in Ramesh et al. (2013)
# determines the strength of bound vortices
# determines the x, z components of the bound vortices
function update_bv(surf::TwoDSurf)
    gamma = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = (surf.a0[1]*(1 + cos(surf.theta[ib])))
        for ia = 1:surf.naterm
            gamma[ib] = gamma[ib] + surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
    end



    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x[ib] + surf.bnd_x[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z[ib] + surf.bnd_z[ib-1])/2.
    end
end

function update_bv(surf::TwoDSurf2DOF)
    gamma = zeros(surf.ndiv)
    for ib = 1:surf.ndiv
        gamma[ib] = (surf.a0[1]*(1 + cos(surf.theta[ib])))
        for ia = 1:surf.naterm
            gamma[ib] = gamma[ib] + surf.aterm[ia]*sin(ia*surf.theta[ib])*sin(surf.theta[ib])
        end
        gamma[ib] = gamma[ib]*surf.uref*surf.c
    end



    for ib = 2:surf.ndiv
        surf.bv[ib-1].s = (gamma[ib]+gamma[ib-1])*(surf.theta[2]-surf.theta[1])/2.
        surf.bv[ib-1].x = (surf.bnd_x[ib] + surf.bnd_x[ib-1])/2.
        surf.bv[ib-1].z = (surf.bnd_z[ib] + surf.bnd_z[ib-1])/2.
    end
end

# Function for calculating the wake rollup
function wakeroll(surf::TwoDSurf, curfield::TwoDFlowField, dt)

    nlev = length(curfield.lev)
    ntev = length(curfield.tev)
    nextv = length(curfield.extv)

    #Clean induced velocities
    for i = 1:ntev
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:nlev
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    for i = 1:nextv
        curfield.extv[i].vx = 0
        curfield.extv[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev; curfield.extv])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(ntev + nlev + nextv)
    wtemp = zeros(ntev + nlev + nextv)
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev); map(q -> q.x, curfield.extv)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev); map(q -> q.z, curfield.extv) ])

    for i = 1:ntev
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vx += utemp[i]
        curfield.lev[i-ntev].vz += wtemp[i]
    end
    for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vx += utemp[i]
        curfield.extv[i-ntev-nlev].vz += wtemp[i]
    end

    #Convect free vortices with their induced velocities
    for i = 1:ntev
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:nlev
            curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    for i = 1:nextv
        curfield.extv[i].x += dt*curfield.extv[i].vx
        curfield.extv[i].z += dt*curfield.extv[i].vz
    end

    return curfield
end

function wakeroll(surf::TwoDSurf2DOF, curfield::TwoDFlowField, dt)

    nlev = length(curfield.lev)
    ntev = length(curfield.tev)
    nextv = length(curfield.extv)

    #Clean induced velocities
    for i = 1:ntev
        curfield.tev[i].vx = 0
        curfield.tev[i].vz = 0
    end

    for i = 1:nlev
        curfield.lev[i].vx = 0
        curfield.lev[i].vz = 0
    end

    for i = 1:nextv
        curfield.extv[i].vx = 0
        curfield.extv[i].vz = 0
    end

    #Velocities induced by free vortices on each other
    mutual_ind([curfield.tev; curfield.lev; curfield.extv])

    #Add the influence of velocities induced by bound vortices
    utemp = zeros(ntev + nlev + nextv)
    wtemp = zeros(ntev + nlev + nextv)
    utemp, wtemp = ind_vel(surf.bv, [map(q -> q.x, curfield.tev); map(q -> q.x, curfield.lev); map(q -> q.x, curfield.extv)], [map(q -> q.z, curfield.tev); map(q -> q.z, curfield.lev); map(q -> q.z, curfield.extv) ])

    for i = 1:ntev
        curfield.tev[i].vx += utemp[i]
        curfield.tev[i].vz += wtemp[i]
    end
    for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vx += utemp[i]
        curfield.lev[i-ntev].vz += wtemp[i]
    end
    for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vx += utemp[i]
        curfield.extv[i-ntev-nlev].vz += wtemp[i]
    end

    #Add influence of gusts on vorticites (only for StepGustDef)

    if (typeof(curfield.velX) == StepGustDef)
       for i = 1:ntev
        curfield.tev[i].vx += curfield.velX(t)
       end
       for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vx += curfield.velX(t)
       end
       for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vx += curfield.velX(t)
       end
    end
    if (typeof(curfield.velZ) == StepGustDef)
       for i = 1:ntev
        curfield.tev[i].vz += curfield.velZ(t)
       end
       for i = ntev+1:ntev+nlev
        curfield.lev[i-ntev].vz += curfield.velZ(t)
       end
       for i = ntev+nlev+1:ntev+nlev+nextv
        curfield.extv[i-ntev-nlev].vz += curfield.velZ(t)
       end
    end

    #Convect free vortices with their induced velocities
    for i = 1:ntev
        curfield.tev[i].x += dt*curfield.tev[i].vx
        curfield.tev[i].z += dt*curfield.tev[i].vz
    end
    for i = 1:nlev
            curfield.lev[i].x += dt*curfield.lev[i].vx
        curfield.lev[i].z += dt*curfield.lev[i].vz
    end
    for i = 1:nextv
        curfield.extv[i].x += dt*curfield.extv[i].vx
        curfield.extv[i].z += dt*curfield.extv[i].vz
    end

    return curfield
end

# Places a trailing edge vortex
function place_tev(surf::TwoDSurf,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc = surf.bnd_z[surf.ndiv]
        else
        xloc = surf.bnd_x[surf.ndiv]+(1./3.)*(field.tev[ntev].x - surf.bnd_x[surf.ndiv])

        zloc = surf.bnd_z[surf.ndiv]+(1./3.)*(field.tev[ntev].z - surf.bnd_z[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end

function place_tev(surf::TwoDSurf2DOF,field::TwoDFlowField,dt)
    ntev = length(field.tev)
    if ntev == 0
        xloc = surf.bnd_x[surf.ndiv] + 0.5*surf.kinem.u*dt
        zloc = surf.bnd_z[surf.ndiv]
        else
        xloc = surf.bnd_x[surf.ndiv]+(1./3.)*(field.tev[ntev].x - surf.bnd_x[surf.ndiv])

        zloc = surf.bnd_z[surf.ndiv]+(1./3.)*(field.tev[ntev].z - surf.bnd_z[surf.ndiv])
    end
    push!(field.tev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))
    return field
end

# Places a leading edge vortex
function place_lev(surf::TwoDSurf,field::TwoDFlowField,dt)
    nlev = length(field.lev)

    le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]

    if (surf.levflag[1] == 0)
        xloc = surf.bnd_x[1] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[1] + 0.5*le_vel_z*dt
    else
        xloc = surf.bnd_x[1]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[1])
        zloc = surf.bnd_z[1]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[1])
    end

    push!(field.lev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))

    return field
end

function place_lev(surf::TwoDSurf2DOF,field::TwoDFlowField,dt)
    nlev = length(field.lev)

    le_vel_x = surf.kinem.u - surf.kinem.alphadot*sin(surf.kinem.alpha)*surf.pvt*surf.c + surf.uind[1]
    le_vel_z = -surf.kinem.alphadot*cos(surf.kinem.alpha)*surf.pvt*surf.c- surf.kinem.hdot + surf.wind[1]

    if (surf.levflag[1] == 0)
        xloc = surf.bnd_x[1] + 0.5*le_vel_x*dt
        zloc = surf.bnd_z[1] + 0.5*le_vel_z*dt
    else
        xloc = surf.bnd_x[1]+(1./3.)*(field.lev[nlev].x - surf.bnd_x[1])
        zloc = surf.bnd_z[1]+(1./3.)*(field.lev[nlev].z - surf.bnd_z[1])
    end

    push!(field.lev,TwoDVort(xloc,zloc,0.,0.02*surf.c,0.,0.))

    return field
end

# Function for updating the positions of the bound vortices
function update_boundpos(surf::TwoDSurf, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x[i] = surf.bnd_x[i] + dt*((surf.pvt*surf.c - surf.x[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + surf.cam[i]*cos(surf.kinem.alpha)*surf.kinem.alphadot)
        surf.bnd_z[i] = surf.bnd_z[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - surf.cam[i]*sin(surf.kinem.alpha)*surf.kinem.alphadot)
    end
    return surf
end

function update_boundpos(surf::TwoDSurf2DOF, dt::Float64)
    for i = 1:surf.ndiv
        surf.bnd_x[i] = surf.bnd_x[i] + dt*((surf.pvt*surf.c - surf.x[i])*sin(surf.kinem.alpha)*surf.kinem.alphadot - surf.kinem.u + surf.cam[i]*cos(surf.kinem.alpha)*surf.kinem.alphadot)
        surf.bnd_z[i] = surf.bnd_z[i] + dt*(surf.kinem.hdot + (surf.pvt*surf.c - surf.x[i])*cos(surf.kinem.alpha)*surf.kinem.alphadot - surf.cam[i]*sin(surf.kinem.alpha)*surf.kinem.alphadot)
    end
    return surf
end

# Function to calculate induced velocities by a set of votices at a target location
function ind_vel(src::Vector{TwoDVort},t_x,t_z)
    #'s' stands for source and 't' for target
    uind = zeros(length(t_x))
    wind = zeros(length(t_x))

    for itr = 1:length(t_x)
        for isr = 1:length(src)
            xdist = src[isr].x - t_x[itr]
            zdist = src[isr].z - t_z[itr]
            distsq = xdist*xdist + zdist*zdist
            uind[itr] = uind[itr] - src[isr].s*zdist/(2*pi*sqrt(src[isr].vc*src[isr].vc*src[isr].vc*src[isr].vc + distsq*distsq))
            wind[itr] = wind[itr] + src[isr].s*xdist/(2*pi*sqrt(src[isr].vc*src[isr].vc*src[isr].vc*src[isr].vc + distsq*distsq))
        end
    end
    return uind, wind
end

# Function determining the effects of interacting vorticies - velocities induced on each other - classical n-body problem
function mutual_ind(vorts::Vector{TwoDVort})
    for i = 1:length(vorts)
        for j = i+1:length(vorts)
            dx = vorts[i].x - vorts[j].x
            dz = vorts[i].z - vorts[j].z
            #source- tar
            dsq = dx*dx + dz*dz

            magitr = 1./(2*pi*sqrt(vorts[j].vc*vorts[j].vc*vorts[j].vc*vorts[j].vc + dsq*dsq))
            magjtr = 1./(2*pi*sqrt(vorts[i].vc*vorts[i].vc*vorts[i].vc*vorts[i].vc + dsq*dsq))

            vorts[j].vx -= dz * vorts[i].s * magjtr
            vorts[j].vz += dx * vorts[i].s * magjtr

            vorts[i].vx += dz * vorts[j].s * magitr
            vorts[i].vz -= dx * vorts[j].s * magitr
        end
    end
    return vorts
end

"""
    controlVortCount(delvort, surf, curfield)

Performs merging or deletion operation on free vortices in order to
control computational cost according to algorithm specified.

Algorithms for parameter delvort

 - delNone

    Does nothing, no vortex count control.

 - delSpalart(limit=500, dist=10, tol=1e-6)

    Merges vortices according to algorithm given in Spalart,
    P. R. (1988). Vortex methods for separated flows.

     - limit: min number of vortices present for merging to occur

     - dist: small values encourage mergin near airfoil, large values in
         wake (see paper)

     - tol: tolerance for merging, merging is less likely to occur for low
        values (see paper)

There is no universal set of parameters that work for all problem. If
using vortex control, test and calibrate choice of parameters

"""
function controlVortCount(delvort :: delVortDef, surf :: TwoDSurf, curfield :: TwoDFlowField)

    if typeof(delvort) == delNone

    elseif typeof(delvort) == delSpalart
        D0 = delvort.dist*surf.c
        V0 = delvort.tol
        if length(curfield.tev) > delvort.limit
            #Check possibility of merging the last vortex with the closest 20 vortices
            gamma_j = curfield.tev[1].s
            d_j = sqrt((curfield.tev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2)
            z_j = sqrt(curfield.tev[1].x^2 + curfield.tev[1].z^2)


            for i = 2:20
                gamma_k = curfield.tev[i].s
                d_k = sqrt((curfield.tev[i].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[i].z- surf.bnd_z[div(surf.ndiv,2)])^2)
                z_k = sqrt(curfield.tev[i].x^2 + curfield.tev[i].z^2)

                fact = abs(gamma_j*gamma_k)*abs(z_j - z_k)/(abs(gamma_j + gamma_k)*(D0 + d_j)^1.5*(D0 + d_k)^1.5)

                if fact < V0
                    #Merge the 2 vortices into the one at k
                    curfield.tev[i].x = (abs(gamma_j)*curfield.tev[1].x + abs(gamma_k)*curfield.tev[i].x)/(abs(gamma_j + gamma_k))
                    curfield.tev[i].z = (abs(gamma_j)*curfield.tev[1].z + abs(gamma_k)*curfield.tev[i].z)/(abs(gamma_j + gamma_k))
                    curfield.tev[i].s += curfield.tev[1].s

                    for j = 1:length(curfield.tev)-1
                        curfield.tev[j] = curfield.tev[j+1]
                    end

                    pop!(curfield.tev)

                    break
                end
            end

            if length(curfield.lev) > delvort.limit
                #Check possibility of merging the last vortex with the closest 20 vortices

                gamma_j = curfield.lev[1].s
                d_j = sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2)
                z_j = sqrt(curfield.lev[1].x^2 + curfield.lev[1].z^2)

                for i = 2:20
                    gamma_k = curfield.lev[i].s
                    d_k = sqrt((curfield.lev[i].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[i].z- surf.bnd_z[div(surf.ndiv,2)])^2)
                    z_k = sqrt(curfield.lev[i].x^2 + curfield.lev[i].z^2)

                    fact = abs(gamma_j*gamma_k)*abs(z_j - z_k)/(abs(gamma_j + gamma_k)*(D0 + d_j)^1.5*(D0 + d_k)^1.5)

                    if fact < V0
                        #Merge the 2 vortices into the one at k
                         curfield.lev[i].x = (abs(gamma_j)*curfield.lev[1].x + abs(gamma_k)*curfield.lev[i].x)/(abs(gamma_j + gamma_k))
                        curfield.lev[i].z = (abs(gamma_j)*curfield.lev[1].z + abs(gamma_k)*curfield.lev[i].z)/(abs(gamma_j + gamma_k))
                        curfield.lev[i].s += curfield.lev[1].s

                        for j = 1:length(curfield.lev)-1
                            curfield.lev[j] = curfield.lev[j+1]
                        end

                        pop!(curfield.lev)

                        break
                    end
                end
            end
        end
    end
end

function controlVortCount(delvort :: delVortDef, surf :: TwoDSurf2DOF, curfield :: TwoDFlowField)

    if typeof(delvort) == delNone

    elseif typeof(delvort) == delSpalart
        D0 = delvort.dist*surf.c
        V0 = delvort.tol
        if length(curfield.tev) > delvort.limit
            #Check possibility of merging the last vortex with the closest 20 vortices
            gamma_j = curfield.tev[1].s
            d_j = sqrt((curfield.tev[1].x - surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[1].z - surf.bnd_z[div(surf.ndiv,2)])^2)
            z_j = sqrt(curfield.tev[1].x^2 + curfield.tev[1].z^2)

            for i = 2:20
                gamma_k = curfield.tev[i].s
                d_k = sqrt((curfield.tev[i].x - surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.tev[i].z - surf.bnd_z[div(surf.ndiv,2)])^2)
                z_k = sqrt(curfield.tev[i].x^2 + curfield.tev[i].z^2)

                fact = abs(gamma_j*gamma_k)*abs(z_j - z_k)/(abs(gamma_j + gamma_k)*(D0 + d_j)^1.5*(D0 + d_k)^1.5)

                if fact < V0
                    #Merge the 2 vortices into the one at k
                    curfield.tev[i].x = (abs(gamma_j)*curfield.tev[1].x + abs(gamma_k)*curfield.tev[i].x)/(abs(gamma_j) + abs(gamma_k))
                    curfield.tev[i].z = (abs(gamma_j)*curfield.tev[1].z + abs(gamma_k)*curfield.tev[i].z)/(abs(gamma_j) + abs(gamma_k))
                    curfield.tev[i].s += curfield.tev[1].s

                    for j = 1:length(curfield.tev)-1
                        curfield.tev[j] = curfield.tev[j+1]
                    end

                    pop!(curfield.tev)

                    break
                end
            end

            if length(curfield.lev) > delvort.limit
                #Check possibility of merging the last vortex with the closest 20 vortices

                gamma_j = curfield.lev[1].s
                d_j = sqrt((curfield.lev[1].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[1].z- surf.bnd_z[div(surf.ndiv,2)])^2)
                z_j = sqrt(curfield.lev[1].x^2 + curfield.lev[1].z^2)

                for i = 2:20
                    gamma_k = curfield.lev[i].s
                    d_k = sqrt((curfield.lev[i].x- surf.bnd_x[div(surf.ndiv,2)])^2 + (curfield.lev[i].z- surf.bnd_z[div(surf.ndiv,2)])^2)
                    z_k = sqrt(curfield.lev[i].x^2 + curfield.lev[i].z^2)

                    fact = abs(gamma_j*gamma_k)*abs(z_j - z_k)/(abs(gamma_j + gamma_k)*(D0 + d_j)^1.5*(D0 + d_k)^1.5)

                    if fact < V0
                        #Merge the 2 vortices into the one at k
                         curfield.lev[i].x = (abs(gamma_j)*curfield.lev[1].x + abs(gamma_k)*curfield.lev[i].x)/(abs(gamma_j) + abs(gamma_k))
                        curfield.lev[i].z = (abs(gamma_j)*curfield.lev[1].z + abs(gamma_k)*curfield.lev[i].z)/(abs(gamma_j) + abs(gamma_k))
                        curfield.lev[i].s += curfield.lev[1].s

                        for j = 1:length(curfield.lev)-1
                            curfield.lev[j] = curfield.lev[j+1]
                        end

                        pop!(curfield.lev)

                        break
                    end
                end
            end
        end
    end
end

function update_kinem(surf::TwoDSurf2DOF, dt, cl, cm)
    #Update previous terms
    surf.kinem.alpha_pr = surf.kinem.alpha
    surf.kinem.h_pr = surf.kinem.h

    surf.kinem.alphadot_pr3 = surf.kinem.alphadot_pr2
    surf.kinem.alphadot_pr2 = surf.kinem.alphadot_pr
    surf.kinem.alphadot_pr = surf.kinem.alphadot

    surf.kinem.hdot_pr3 = surf.kinem.hdot_pr2
    surf.kinem.hdot_pr2 = surf.kinem.hdot_pr
    surf.kinem.hdot_pr = surf.kinem.hdot

    surf.kinem.alphaddot_pr3 = surf.kinem.alphaddot_pr2
    surf.kinem.alphaddot_pr2 = surf.kinem.alphaddot_pr

    surf.kinem.hddot_pr3 = surf.kinem.hddot_pr2
    surf.kinem.hddot_pr2 = surf.kinem.hddot_pr

    # Calculate hddot and alphaddot from forces based on 2DOF structural system
    m11 = 2./surf.c
    m12 = -surf.strpar.x_alpha*cos(surf.kinem.alpha)
    m21 = -2.*surf.strpar.x_alpha*cos(surf.kinem.alpha)/surf.c
    m22 = surf.strpar.r_alpha*surf.strpar.r_alpha

    R1 = 4*surf.strpar.kappa*surf.uref*surf.uref*cl/(pi*surf.c*surf.c) - 2*surf.strpar.w_h*surf.strpar.w_h*(surf.strpar.cubic_h_1*surf.kinem.h + surf.strpar.cubic_h_3*surf.kinem.h^3)/surf.c - surf.strpar.x_alpha*sin(surf.kinem.alpha)*surf.kinem.alphadot*surf.kinem.alphadot

    R2 = 8*surf.strpar.kappa*surf.uref*surf.uref*cm/(pi*surf.c*surf.c) - surf.strpar.w_alpha*surf.strpar.w_alpha*surf.strpar.r_alpha*surf.strpar.r_alpha*(surf.strpar.cubic_alpha_1*surf.kinem.alpha + surf.strpar.cubic_alpha_3*surf.kinem.alpha^3)

    surf.kinem.hddot_pr = (1/(m11*m22-m21*m12))*(m22*R1-m12*R2)
    surf.kinem.alphaddot_pr = (1/(m11*m22-m21*m12))*(-m21*R1+m11*R2)

    #Kinematics are updated according to the 2DOF response
    surf.kinem.alphadot = surf.kinem.alphadot_pr + (dt/12.)*(23*surf.kinem.alphaddot_pr - 16*surf.kinem.alphaddot_pr2 + 5*surf.kinem.alphaddot_pr3)
    surf.kinem.hdot = surf.kinem.hdot_pr + (dt/12)*(23*surf.kinem.hddot_pr-16*surf.kinem.hddot_pr2 + 5*surf.kinem.hddot_pr3)

    surf.kinem.alpha = surf.kinem.alpha_pr + (dt/12)*(23*surf.kinem.alphadot_pr-16*surf.kinem.alphadot_pr2 + 5*surf.kinem.alphadot_pr3)
    surf.kinem.h = surf.kinem.h_pr + (dt/12)*(23*surf.kinem.hdot_pr - 16*surf.kinem.hdot_pr2 + 5*surf.kinem.hdot_pr3)

    return surf
end
