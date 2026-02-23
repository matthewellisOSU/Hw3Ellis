from math import sqrt, pi

def gamma_rec(alpha):
    """
    Computes gamma for alpha that is an integer or half integer
    base cases:
    gamma(1) = 1
    gamma(1/2) = sqrt(pi)

    Reoccurance:
    gamma(a) = (a-1)*gamma(a-1)
    :param alpha: value where gamma is evaluated
    :return: gamma
    """
    if abs(alpha - 1.0) < 1e-12:
        return 1.0
    if abs(alpha - 0.5) < 1e-12:
        return sqrt(pi)
    return (alpha - 1.0) * gamma_rec(alpha - 1.0)

def km_constant(m):
    """
    computes km constant
    :param m: int degrees of freedom
    :return: km constant
    """
    return gamma_rec((m + 1) / 2.0) / (sqrt(m * pi) * gamma_rec(m / 2.0))

def t_integrand(u, m):
    """
    computes t_integrand
    :param u: integration variable
    :param m: int degrees of freedom
    :return: integrand value
    """
    return (1.0 + (u * u) / m) ** (-(m + 1) / 2.0)

def simpson(fn, a, b, m, N=600):
    """
    Uses Simpson's 1/3 rule to integrate the normal pdf between limits a and b.
    :param fn: callback function
    :param a: left limit
    :param b: right limit
    :param m: int degrees of freedom
    :param N: # subintervals
    :return: approximate integration value
    """

    if N % 2 == 1:
        N += 1

    h = (b - a) / N
    total = 0.0

    for i in range(N + 1):
        u = a + i * h
        fu = fn(u, m)

        if i == 0 or i == N:
            total += fu
        elif i % 2 == 1:
            total += 4.0 * fu
        else:
            total += 2.0 * fu

    return (h / 3.0) * total

def t_cdf(z, m):
    """
    Computes CDF of F(z) w/ m degrees of freedom
    :param z: value to compute CDF
    :param m: int degrees of freedom
    :return: P
    """
    Km = km_constant(m)

    if z == 0.0:
        return 0.5

    if z > 0:
        area = simpson(t_integrand, 0.0, z, m, N=600)
        return 0.5 + Km * area
    else:
        area = simpson(t_integrand, 0.0, -z, m, N=600)
        return 0.5 - Km * area

def main():
    """
    Main program that asks user for m and z and F(z) values
    :return: 
    """
    print("HW3b: t-distribution CDF F(z)")
    m = int(input("Enter degrees of freedom m: "))
    z = float(input("Enter z value: "))

    p = t_cdf(z, m)
    print("F({:0.4f}) for m={} = {:0.6f}".format(z, m, p))


if __name__ == "__main__":
    main()

#ChatGPT was consulted for logic and correctness