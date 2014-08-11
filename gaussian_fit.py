

def gaussian_no_background(x, area, mu, sigma):
    from pylab import sqrt, pi, exp
    return area/(sigma*sqrt(2*pi))*exp(-((x-mu)/sigma)**2/2)
  
def gaussian_linear_background(x, area, mu, sigma, a, b):
    from pylab import sqrt, pi, exp
    return area/(sigma*sqrt(2*pi))*exp(-((x-mu)/sigma)**2/2) + a*x + b
