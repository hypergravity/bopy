# -*- coding: utf-8 -*-
"""
@author: cham
Created on Sat Jul  4 21:05:41 2015
"""
#%%
from numpy import *
Nobs = 20
x_true = random.uniform(0,10, size=Nobs)    # size(x_true): 20
y_true = random.uniform(-1,1, size=Nobs)    # size(y_true): 20
alpha_true = 0.5                            # alpha_true:   0.5
beta_x_true = 1.0                           # beta_x_true:  1.0
beta_y_true = 10.0                          # beta_y_true:  10.0
eps_true = 0.5                              # eps_true:     0.5
z_true = alpha_true + beta_x_true*x_true + beta_y_true*y_true   # model
z_obs = z_true + random.normal(0, eps_true, size=Nobs)          # std = eps_true

#%% x|y vs z_true
from matplotlib import pyplot as plt
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.scatter(x_true, z_true, c=y_true, marker='o')    # x_true vs z_obs
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Z')
plt.subplot(1,2,2)
plt.scatter(y_true, z_true, c=x_true, marker='o')    # y_true vs z_obs
plt.colorbar()
plt.xlabel('Y')
plt.ylabel('Z')

#%% x|y vs z_obs
#%matplotlib inline
from matplotlib import pyplot as plt
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.scatter(x_true, z_obs, c=y_true, marker='o')    # x_true vs z_obs
plt.colorbar()
plt.xlabel('X')
plt.ylabel('Z')
plt.subplot(1,2,2)
plt.scatter(y_true, z_obs, c=x_true, marker='o')    # y_true vs z_obs
plt.colorbar()
plt.xlabel('Y')
plt.ylabel('Z')

#%%
def lnprior(p):
    # The parameters are stored as a vector of values, so unpack them
    alpha,betax,betay,eps = p
    # We're using only uniform priors, and only eps has a lower bound
    if eps <= 0:
        return -inf
    return 0

def lnlike(p, x, y, z):
    alpha,betax,betay,eps = p
    model = alpha + betax*x + betay*y
    # the likelihood is sum of the lot of normal distributions
    denom = power(eps,2)
    lp = -0.5*sum(power((z - model),2)/denom + log(denom) + log(2*pi))
    return lp

def lnprob(p, x, y, z):
    lp = lnprior(p)
    if not isfinite(lp):
        return -inf
    return lp + lnlike(p, x, y, z)

#%%
import scipy.optimize as opt
nll = lambda *args: -lnlike(*args)
result = opt.minimize(nll, [alpha_true, beta_x_true, beta_y_true, eps_true],args=(x_true, y_true, z_obs))
print result['x']


#%%
Nwalker,Ndim = 50,4
p0 = [result['x']+1.e-4*random.randn(Ndim) for i in range(Nwalker)]

#%%
from multiprocessing.dummy import Pool as ThreadPool
pool_inst = ThreadPool(4)
import emcee
#rc = parallel.client()
sampler = emcee.EnsembleSampler(Nwalker,Ndim,lnprob,args=(x_true,y_true,z_obs),pool=pool_inst)
pos,prob,state = sampler.run_mcmc(p0, 500)

#%%
#%pylab inline
res=plot(sampler.chain[:,:,0].T, '-', color='k', alpha=0.3)
axhline(alpha_true, color='blue')

#%%
sampler.reset()
pos,prob,state = sampler.run_mcmc(pos, 1000)

#%%
m_alpha,m_betax,m_betay,m_eps = median(sampler.flatchain, axis=0)

plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.plot(x_true, z_obs-m_alpha-m_betay*y_true, 'o')
plt.xlabel('X')
plt.ylabel('Z - alpha - beta_y y')
# Now plot the model
xx = array([x_true.min(), x_true.max()])
plt.plot(xx, xx*m_betax)
plt.plot(xx, xx*m_betax + m_eps, '--', color='k')
plt.plot(xx, xx*m_betax - m_eps, '--', color='k')
plt.subplot(1,2,2)
plt.plot(y_true, z_obs-m_alpha-m_betax*x_true, 'o')
plt.xlabel('Y')
plt.ylabel('Z - alpha - beta_x x')
yy = array([y_true.min(), y_true.max()])
plt.plot(yy, yy*m_betay)
plt.plot(yy, yy*m_betay + m_eps, '--', color='k')
plt.plot(yy, yy*m_betay - m_eps, '--', color='k')

#%% triangle_plot
import triangle
triangle.corner(sampler.flatchain)
tmp = triangle.corner(sampler.flatchain, labels=['alpha','betax','betay','eps'],truths=[alpha_true, beta_x_true, beta_y_true, eps_true])

