def quad_darkening(x, a)
    return a[0] * (1.0 - a[1] - a[2] + a[1] * x + a[2] * x^2)
