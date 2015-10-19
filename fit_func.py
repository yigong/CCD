from numpy import sqrt, pi, exp
import numpy as np
from scipy.optimize import curve_fit

def myGaussian(x, mu, sigma):
    return exp(-((x-mu)/sigma)**2/2)

def Gaussian(x, a0, mu, sigma):
    return a0*exp(-((x-mu)/sigma)**2/2)

def Gaussian_area(x, area, mu, sigma):
    return area/(sigma*sqrt(2*pi))*exp(-((x-mu)/sigma)**2/2)

def Gaussian_const_bg(x, area, mu, sigma, b):
    return area/(sigma*sqrt(2*pi))*exp(-((x-mu)/sigma)**2/2)  + b

def Gaussian_double(x, area1, mu1, sigma1, area2, mu2, sigma2):
    return area1/(sigma1*sqrt(2*pi))*exp(-((x-mu1)/sigma1)**2/2) + \
           area2/(sigma2*sqrt(2*pi))*exp(-((x-mu2)/sigma2)**2/2)

def fit_Gaussian(x, y, a0, mu, sigma):
    val, covMX = curve_fit(Gaussian, x, y, p0=[a0, mu, sigma])
    cov = np.diag(covMX)
    return val, cov

def fit_Gaussian_area(x, y, area, mu, sigma):
    val, covMX = curve_fit(Gaussian_area, x, y, p0=[area, mu, sigma])
    cov = np.diag(covMX)
    return val, cov


def Cauchy(x, a0, x0, g):
    return a0/(np.pi * g * (1+((x-x0)/g)**2))

def fit_Cauchy(x, y, a0, x0, g):
    val, covMX = curve_fit(Cauchy, x, y, p0=[a0, x0, g])
    cov = np.diag(covMX)
    return val, cov

def myCauchy(x, x0, g):
    return 1/(1 + (((x-x0)/g)**2))

def fit_myCauchy(x, y, x0, g):
    val, covMX = curve_fit(myCauchy, x, y, p0=[x0, g])
    cov = np.diag(covMX)
    return val, cov

