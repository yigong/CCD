from numpy import sqrt, pi, exp

def gauss(x, a0, mu, sigma):
    return a0*exp(-((x-mu)/sigma)**2/2)

def gauss_area(x, area, mu, sigma):
    return area/(sigma*sqrt(2*pi))*exp(-((x-mu)/sigma)**2/2)

def gauss_const_bg(x, area, mu, sigma, b):
    return area/(sigma*sqrt(2*pi))*exp(-((x-mu)/sigma)**2/2)  + b

def doubleGaussianFit(x, area1, mu1, sigma1, area2, mu2, sigma2):
    return area1/(sigma1*sqrt(2*pi))*exp(-((x-mu1)/sigma1)**2/2) + \
           area2/(sigma2*sqrt(2*pi))*exp(-((x-mu2)/sigma2)**2/2)
