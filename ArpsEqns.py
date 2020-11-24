import numpy as np

def SecantToTangent(di_s, b):

    if b == 0:
        di_t = di_s
    else:
        di_t = 1-np.exp(-((1-di_s)**(-b)-1)/b)

    return di_t

def SecantToNominal(di_s, b):

    if b == 0:
        ai = -np.log(1-di_s)
    elif b == 1:
        ai = di_s / (1-di_s)
    else:
        ai = (1/b)*((1-di_s)**(-b)-1)

    return ai

def NominalToSecant(ai, b):

    if b == 0:
        di_s = 1-np.exp(-ai)
    elif b == 1:
        di_s = ai/(1+ai)
    else:
        di_s = 1-(1+b*ai)**(-1/b)

    return di_s

def TangentToNominal(di_t, b):

    ai = -np.log(1-di_t)

    return ai

def TangentToSecant(di_t, b):

    di_s = NominalToSecant(TangentToNominal(di_t, b),b)

    return di_s

def NominalToTangent(ai, b):

    di_t = SecantToTangent(NominalToSecant(ai,b),b)

    return di_t

def ArpsRate(qi, di_s, b, t, dmin=0, bmin=0):
    """output: rate at time t in units of vol/d
       qi - initial rate (vol/d)
       di_s - initial decline (%/yr - secant)
       b - b-factor
       t - time (days)
       dmin - optional modified arps minimum decline rate (%/yr - tangent)
       bmin - optional modified arps b-factor"""

    def RateEqn(q,a,b,t):
        if b == 0:
            r = q * np.exp(-a*t)
        elif b == 1:
            r = q / (1+a*t)
        else:
            r = q/((1+b*a*t)**(1/b))
        return r

    ai = SecantToNominal(di_s, b)/365

    if not (dmin == 0 or b == 0):
        # Get nominal terminal decline, but note that dmin is in terms of tangential decline
        aimin = TangentToNominal(dmin, bmin)/365
        qmin = qi * (aimin/ai)**(1/b)
        tmin = ((qi/qmin)**b-1)/(b*ai)
    else:
        tmin = -1

    if t < tmin or tmin == -1:
        rate = RateEqn(qi, ai, b, t)
    else:
        rate = RateEqn(qmin, aimin, bmin, (t-tmin))

    return rate

def ArpsNp(qi, di_s, b, t, dmin=0, bmin=0):
    """output: cumulative volume from time zero to time t
       qi - initial rate (vol/d)
       di_s - initial decline (%/yr - secant)
       b - b-factor
       t - time (days)
       dmin - optional modified arps minimum decline rate (%/yr - tangent)
       bmin - optional modified arps b-factor"""

    def RateEqn(q,a,b,t):
        if b == 0:
            r = q * np.exp(-a*t)
        elif b == 1:
            r = q / (1+a*t)
        else:
            r = q/((1+b*a*t)**(1/b))
        return r

    def CumEqn(q,a,b,qt):
        if b == 0:
            cum = (q-qt)/a
        elif b == 1:
            cum = (q/a)*np.log(q/qt)
        else:
            cum = (q**b)/(a*(1-b))*(q**(1-b)-qt**(1-b))
        return cum

    ai = SecantToNominal(di_s, b)/365

    if not (dmin == 0 or b == 0):
        # Get nominal terminal decline, but note that dmin is in terms of tangential decline
        aimin = TangentToNominal(dmin, bmin)/365
        qmin = qi * (aimin/ai)**(1/b)
        tmin = ((qi/qmin)**b-1)/(b*ai)
        npi = CumEqn(qi, ai, b, qmin)
    else:
        tmin = -1

    if t < tmin or tmin == -1:
        rate = RateEqn(qi, ai, b, t)
        Q = CumEqn(qi, ai, b, rate)
    else:
        rate = RateEqn(qmin, aimin, bmin, (t-tmin))
        Q = npi + CumEqn(qmin, aimin, bmin, rate)

    return Q

def ArpsLife(qi, di_s, b, qab, dmin=0, bmin=0, max_life = 50*365):
    """output: time to abandonment rate in days
       qi - initial rate (vol/d)
       di_s - initial decline (%/yr - secant)
       b - b-factor
       qab - abandonment rate (vol/d)
       dmin - optional modified arps minimum decline rate (%/yr - tangent)
       bmin - optional modified arps b-factor
       max_life - optional life cutoff"""

    def LifeEqn(q, a, b, qab):

        if b == 0:
            life = np.log(q/qab)/a
        else:
            life = ((q/qab)**b-1)/(b*a)

        return life

    ai = SecantToNominal(di_s, b) / 365

    if not (dmin == 0 or b == 0):
        # Get nominal terminal decline, but note that dmin is in terms of tangential decline
        aimin = TangentToNominal(dmin, bmin) / 365
        qmin = qi * (aimin / ai) ** (1 / b)
        tmin = ((qi / qmin) ** b - 1) / (b * ai)
    else:
        qmin = -1
        tmin = -1

    if qab > qmin or tmin == -1:
        life = LifeEqn(qi, ai, b, qab)
    else:
        life = tmin + LifeEqn(qmin, aimin, bmin, qab)

    if life > max_life:
        life = max_life

    return life

def ArpsEUR(qi, di_s, b, qab, dmin=0, bmin=0, max_life = 50*365):
    """output: EUR
       qi - initial rate (vol/d)
       di_s - initial decline (%/yr - secant)
       b - b-factor
       qab - abandonment rate (vol/d)
       dmin - optional modified arps minimum decline rate (%/yr - tangent)
       bmin - optional modified arps b-factor
       max_life - optional life cutoff"""

    def LifeEqn(q, a, b, qab):

        if b == 0:
            life = np.log(q/qab)/a
        else:
            life = ((q/qab)**b-1)/(b*a)

        return life

    ai = SecantToNominal(di_s, b) / 365

    if not (dmin == 0 or b == 0):
        # Get nominal terminal decline, but note that dmin is in terms of tangential decline
        aimin = TangentToNominal(dmin, bmin) / 365
        qmin = qi * (aimin / ai) ** (1 / b)
        tmin = ((qi / qmin) ** b - 1) / (b * ai)
    else:
        qmin = -1
        tmin = -1

    if qab > qmin or tmin == -1:
        life = LifeEqn(qi, ai, b, qab)
        if life > max_life:
            life = max_life
        eur = ArpsNp(qi, di_s, b, life, dmin, bmin)
    else:
        life = tmin + LifeEqn(qmin, aimin, bmin, qab)
        if life > max_life:
            life = max_life
        eur = ArpsNp(qi, di_s, b, life, dmin, bmin)

    return eur

def ArpsBDEG(qi, di_s, b, qab, target_percent, dmin=0, bmin=0, max_life = 50*365, ramp_cum = 0):
    """output: b-factor that results in target percent of original eur by holding all inputs except b constant
       NOTE: nominal/tangential decline is held constant, therefore as b changes, it also changes secant decline
       qi - initial rate (vol/d)
       di_s - initial decline (%/yr - secant)
       b - b-factor
       qab - abandonment rate (vol/d)
       target_percent - percent of eur of initial inputs
       dmin - optional modified arps minimum decline rate (%/yr - tangent)
       bmin - optional modified arps b-factor
       max_life - optional life cutoff
       ramp_cum - optional total cumulative volume contained in the ramp phase"""

    original_eur = ArpsEUR(qi, di_s, b, qab, dmin, bmin, max_life) + ramp_cum
    target_eur = original_eur * target_percent
    original_di_t = SecantToTangent(di_s, b)

    # Initial guess and direction
    errorcount = 0
    bdeg = b
    incr = -0.1
    temp_di_s = TangentToSecant(original_di_t, bdeg)
    temp_eur = ArpsEUR(qi, temp_di_s, bdeg, qab, dmin, bmin, max_life) + ramp_cum

    while round(target_eur, 1) != round(temp_eur, 1):

        bdeg = bdeg + incr

        if bdeg < 0:
            bdeg = 0
            errorcount += 1
        if errorcount == 3:
            return bdeg

        temp_di_s = TangentToSecant(original_di_t, bdeg)
        temp_eur = ArpsEUR(qi, temp_di_s, bdeg, qab, dmin, bmin, max_life) + ramp_cum

        if not ((temp_eur > target_eur and incr < 0) or (temp_eur < target_eur and incr > 0)):
            # Getting close, change direction of guess and make smaller
            incr = incr * -0.1

    return bdeg


