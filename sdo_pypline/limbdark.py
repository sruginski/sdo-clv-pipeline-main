def quad_darkening(x, a, b, c):
    return a * (1.0 - b * (1.0 - x) - c * (1.0 - x)**2)

# def quad_darkening(x, a):
#     return a[0] * (1.0 - a[1] * (1.0 - x) - a[2] * (1.0 - x)**2)
