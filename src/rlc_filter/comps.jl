
function TimeDependentInductor(; name, i_start=0.0)
    @named oneport = OnePort(; i_start=i_start)
    @unpack v, i = oneport
    sts = @variables Lt(t)
    eqs = [
        D(i) ~ 1 / Lt * v,
    ]
    extend(ODESystem(eqs, t, sts, []; name=name), oneport)
end

function TimeDependentCapacitor(; name, v_start=0.0)
    @named oneport = OnePort(; v_start=v_start)
    @unpack v, i = oneport
    # pars = @parameters C = C
    sts = @variables Ct(t)
    eqs = [
        D(v) ~ i / Ct,
    ]
    extend(ODESystem(eqs, t, sts, []; name=name), oneport)
end

function TimeDependentResistor(; name)
    @named oneport = OnePort()
    @unpack v, i = oneport
    # pars = @parameters R = R
    sts = @variables Rt(t)
    eqs = [
        v ~ i * Rt,
    ]
    extend(ODESystem(eqs, t, sts, []; name=name), oneport)
end

@named o = Constant(; k=1.0)
@named o = Sine(; frequency=1)
