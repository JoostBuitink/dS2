import numpy as np

def gQ_fun(Q, alpha, beta, gamma):
    '''
    Discharge sensitivy function (sensitivy of discharge to changes in
    storage).

    Input
    -----
    alpha: float
        alpha parameter
    beta: float
        beta parameter
    gamma: float
        gamma parameter
    Q: float
        Discharge

    Output
    ------
    gQ: float
        Result of np.exp(alpha + beta * np.log(Q) + (gamma / Q))
        Is essentially the same as a * Q ** beta * np.exp(gamma / Q), where
        a = np.exp(alpha), but the exponential function is a lot faster
    '''
    gQ = np.exp(alpha + beta * np.log(Q) + (gamma / Q))
    gQ = np.fmax(gQ, 1e-200)
    return gQ

def dQdt_fun(Q, P, ET, alpha, beta, gamma):
    return gQ_fun(Q, alpha, beta, gamma) * (P - ET - Q)

def flexible_RK4(P, ET, Q, alpha, beta, gamma, dt, LB,
                 gQ_mdiff=1.25, dt_reduction=.15, min_dt=5, max_dt=100):
    '''
    Numerical solver based on RK4, but with a flexible time stepping scheme
    implemented. Timestep is reduced based on the gQ value (if gQ > 1) or
    based on the ratio between gQt-1 and gQt.

    Input
    -----
    P: float or numpy.array
        Precipitation at time t in mm h-1
    ET: float or numpy.array
        Evaporation at time t in mm h-1
    Q: float or numpy.array
        Discharge at time t-1 in mm h-1
    alpha: float or numpy.array
        gQ parameter value
    beta: float or numpy.array
        gQ parameter value
    gamma: float or numpy.array
        gQ parmater value
    dt: int
        Timestep in hours
    LB: float
        Lower boundary to prevent negative discharge values (fraction of Q)
    gQ_midff: float
        Maximum allowed relative difference between g(Q_t) and g(Q_t-1)
    dt_reduction: float
        Factor determining the how many extra timesteps are required
        (function of gQ_mdiff)
    min_dt: int
        Minimum extra timesteps if time step needs to be reduced
    max_dt: int
        Maximum extra timesteps if time step needs to be reduced

    Output
    ------
    Qout: float or numpy.array
        Discharge value at time step t+dt

    See Also
    --------
    :func:`normal_RK4`: Runge-Kutta 4 implementation
    :func:`gQ_fun`: discharge sensitivity function

    '''
    Qout, gQ0 = normal_RK4(P, ET, Q, alpha, beta, gamma, dt, LB)
    gQout = gQ_fun(Qout, alpha, beta, gamma)

    gQ_diff = np.max(np.abs(gQout - gQ0) / np.minimum(gQout, gQ0))
    number_dt = 1

    if np.any(gQout >= 1):
        number_dt = int(max(min_dt, min(np.max(gQout) * 10, max_dt)))
        Qprev = Q
        for i in range(number_dt):
            Qprev, tmp = normal_RK4(P, ET, Qprev, alpha, beta, gamma, dt/number_dt, LB)
        Qout = Qprev
    elif gQ_diff > gQ_mdiff:
        number_dt = int(max(min_dt, min(gQ_diff ** dt_reduction, max_dt)))
        Qprev = Q
        for i in range(number_dt):
            Qprev, tmp = normal_RK4(P, ET, Qprev, alpha, beta, gamma, dt/number_dt, LB)
        Qout = Qprev

    return Qout, number_dt

def normal_RK4(P, ET, Q, alpha, beta, gamma, dt, LB):
    '''
    Modified Runge-Kutta 4 scheme, so that the intermediate value is presented
    as the absolute value, in stead of the change. Also returns the gQ value
    of gQ(t-1) '''

    cor = Q * LB
    gQs0 = gQ_fun(Q, alpha, beta, gamma)
    s1 = Q + dt/2 * gQs0 * (P - ET - Q)
    s1 = np.fmax(s1, cor)
    s2 = Q + dt/2 * gQ_fun(s1, alpha, beta, gamma) * (P - ET - s1)
    s2 = np.fmax(s2, cor)
    s3 = Q + dt   * gQ_fun(s2, alpha, beta, gamma) * (P - ET - s2)
    s3 = np.fmax(s3, cor)
    s4 = Q + dt   * gQ_fun(s3, alpha, beta, gamma) * (P - ET - s3)

    Qout = 1/3*s1 + 2/3*s2 + 1/3*s3 + 1/6*s4 - 1/2*Q
    Qout = np.fmax(Qout, cor)
    return Qout, gQs0


def cashkarp_wrapper(P, ET, Q, alpha, beta, gamma, dt=1, LB=1e-10,
               dt_factor = 1000, min_dt = 5, max_dt = 100):
    ''' Wrapper for the Runge-Kutta Cash-Karp solver, which uses the fourth and
    fifth order solutions to estimate the numerical error, and the estimate the
    required fine timestep. '''

    fifth_sol, fourth_sol = rk_cashkarp(P=P, ET=ET, Q=Q, alpha=alpha, beta=beta, gamma=gamma, dt=dt)

    error = np.max(abs((fifth_sol - fourth_sol) / fifth_sol))
    extradt = 1

    if error >= 1e-3 or np.isnan(error):
        extradt = int(min(max(min_dt, error * dt_factor), max_dt)) if not np.isnan(error) else max_dt
        Qprev = Q
        for j in range(extradt):
            Qprev = rk_cashkarp(Q=Qprev, P=P, ET=ET, alpha=alpha, beta=beta, gamma=gamma, dt=dt/extradt)[0]
        fifth_sol = np.fmax(Qprev, 1e-10)

    return fifth_sol, extradt

def rk_cashkarp(P, ET, Q, alpha, beta, gamma, dt=1):
    '''
    Cash-Karp method for solving the differential equation, where both a fourth
    and fifth order solutions are calculated. This difference can be used to
    estimate the numerical error, which in turn can be used to select the a
    finer timestep to reduce this error.
    '''
    Qs1 = Q
    s1 = dQdt_fun(Q=Qs1, P=P, ET=ET, alpha=alpha, beta=beta, gamma=gamma)

    Qs2 = Q + dt * (1/5 * s1)
    s2 = dQdt_fun(Q=Qs2, P=P, ET=ET, alpha=alpha, beta=beta, gamma=gamma)

    Qs3 = Q + dt * (3/40 * s1 + 9/40 * s2)
    s3 = dQdt_fun(Q=Qs3, P=P, ET=ET, alpha=alpha, beta=beta, gamma=gamma)

    Qs4 = Q + dt * (3/10 * s1 - 9/10 * s2 + 6/5 * s3)
    s4 = dQdt_fun(Q=Qs4, P=P, ET=ET, alpha=alpha, beta=beta, gamma=gamma)

    Qs5 = Q + dt * (-11/54 * s1 + 5/2 * s2 - 70/27 * s3 + 35/27 * s4)
    s5 = dQdt_fun(Q=Qs5, P=P, ET=ET, alpha=alpha, beta=beta, gamma=gamma)

    Qs6 = Q + dt * (1631/55296 * s1 + 175/512 * s2 + 575/13824 * s3 + 44275/110592 * s4 + 253/4096 * s5)
    s6 = dQdt_fun(Q=Qs6, P=P, ET=ET, alpha=alpha, beta=beta, gamma=gamma)

    fifth_sol = Q + dt * (37/378 * s1 + 250/621 * s3 + 125/594 * s4 + 512/1771 * s6)
    fourth_sol = Q + dt * (2825/27648 * s1 + 18575/48384 * s3 + 13525/55296 * s4 + 277/14336 * s5 + 1/4 * s6)

    return fifth_sol, fourth_sol

