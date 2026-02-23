from math import sqrt, pi, exp


def gpdf(x, mu, sig):
    """
    Computes the Gaussian probability density function value.

    :param x: value where pdf is evaluated
    :param mu: population mean
    :param sig: population standard deviation
    :return: pdf value at x
    """
    return (1.0 / (sig * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sig) ** 2)


def simpson(mu, sig, a, b, N=200):
    """
    Uses Simpson's 1/3 rule to integrate the normal pdf between limits a and b.

    :param mu: population mean
    :param sig: population standard deviation
    :param a: left hand limit
    :param b: right hand limit
    :param N: number of intervals
    :return: numerical integral result
    """
    if N % 2 == 1:
        N += 1

    h = (b - a) / N
    total = 0.0

    for i in range(N + 1):
        x = a + i * h
        fx = gpdf(x, mu, sig)

        if i == 0 or i == N:
            total += fx
        elif i % 2 == 1:
            total += 4 * fx
        else:
            total += 2 * fx

    return (h / 3) * total


def prob_single(mu, sig, c, GT):
    """
    Computes a single sided probability P(x<c) or P(x>c).

    :param mu: population mean
    :param sig: population standard deviation
    :param c: cutoff value
    :param GT: True if computing P(x>c), False if P(x<c)
    :return: probability value
    """
    a = mu - 5 * sig
    p_less = simpson(mu, sig, a, c, 400)

    if GT:
        return 1 - p_less
    else:
        return p_less


def prob_double(mu, sig, c, inside):
    """
    Computes a double sided probability around the mean.

    :param mu: population mean
    :param sig: population standard deviation
    :param c: boundary value
    :param inside: True for inside band, False for outside band
    :return: probability value
    """
    d = abs(c - mu)
    a = mu - d
    b = mu + d

    p_inside = simpson(mu, sig, a, b, 400)

    if inside:
        return p_inside
    else:
        return 1 - p_inside


def secant(f, x0, x1, tol=1e-6, maxiter=30):
    """
    Uses the Secant method to find a root of a function.

    :param f: function to solve
    :param x0: first guess
    :param x1: second guess
    :param tol: tolerance for stopping
    :param maxiter: maximum iterations
    :return: root estimate and iteration count
    """
    it = 0

    while it < maxiter:
        f0 = f(x0)
        f1 = f(x1)

        if (f1 - f0) == 0:
            break

        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)

        if abs(x2 - x1) < tol:
            return x2, it

        x0 = x1
        x1 = x2
        it += 1

    return x1, it


def main():
    """
    Main program that allows user to compute probabilities or solve for c
    using single or double sided normal distributions.
    """
    mu = float(input("mu: "))
    sig = float(input("sigma: "))

    mode = input("single or double (s/d): ").lower()
    solve = input("specify c or P (c/p): ").lower()

    if mode == "s":
        side = input("> or <: ")
        GT = True if side == ">" else False

        if solve == "c":
            c = float(input("c: "))
            p = prob_single(mu, sig, c, GT)
            print("P =", p)
        else:
            p_target = float(input("P: "))

            def f(cval):
                return prob_single(mu, sig, cval, GT) - p_target

            x0 = mu
            x1 = mu + sig if GT else mu - sig

            c_sol, it = secant(f, x0, x1)
            print("c =", c_sol)
            print("iterations =", it)

    else:
        inside = input("inside band? (y/n): ").lower() == "y"

        if solve == "c":
            c = float(input("c: "))
            p = prob_double(mu, sig, c, inside)
            print("P =", p)
        else:
            p_target = float(input("P: "))

            def f(cval):
                return prob_double(mu, sig, cval, inside) - p_target

            x0 = mu
            x1 = mu + sig

            c_sol, it = secant(f, x0, x1)
            print("c =", c_sol)
            print("iterations =", it)


if __name__ == "__main__":
    main()

#ChatGPT was consulted for logic and correctness