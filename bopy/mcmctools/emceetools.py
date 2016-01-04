# -*- coding: utf-8 -*-
"""
@author: cham
Created on Mon Jan  4 20:08:23 2016
"""

import numpy as np
import emcee
import scipy.optimize as opt


def emcee_general_run(fun_prior,
                      fun_model,
                      fun_prob,
                      data,
                      theta0='scipyopt',
                      theta_scale,
                      param_labels,
                      output_dir,
                      output_prefix,
                      N=(50, 3, 1, 5000, 1000)):
    """ use this to perform a general MCMC/emcee run

    Parameters:
    -----------
    fun_prior: function
        The function returns prior of theta.

    fun_model: function
        The function returns the model results.

    fun_prob: function
        The function returns the likelihood

    data: array or list of arrays
        The list contains data. This should fit the fun_model.

    theta0: {tuple/list/array or 'scipyopt',}
        Initial guess for theta0.
        The default value is 'scipyopt', which uses scipy.opt function to
        generate a best fit theta0.

    theta_scale: array
        To get a reasonable parameter range, theta should be scaled
        by this array.

    output_dir: string
        The output directory path.

    output_prefix: string
        The prefix of output files.

    N: tuple (N_walker, N_dim, N_threads, N_burnin, N_use)
        The set of integer numbers used in the MCMC run.
        The default value is (50, 3, 1, 500, 100), whenever fits the data
        or not.
    """
    # unpack parameters
    N_walker, N_dim, Nthreads, N_burnin, N_use = N
    output_parameter_check = '%s/%s_burnin_%d_use_%d_scale_%s_pcheck.png' % \
        (output_dir, output_prefix, fit_scale, N_burnin, N_use)
    output_corner = '%s/%s_burnin_%d_use_%d_scale_%s_corner.png' % \
        (output_dir, output_prefix, fit_scale, N_burnin, N_use)
    output_chaindata = '%s/%s_burnin_%d_use_%d_scale_%s_chaindata.png' % \
        (output_dir, output_prefix, fit_scale, N_burnin, N_use)

    # initial guess
    if theta0 is 'scipyopt':
        nll = lambda *args: - fun_prob(*args)
        p0 = [1., 1., 1., 1.]
        result = opt.minimize(nll, p0, args=data)
        p0 = result['x']
    else:
        p0 = theta0
    p0 = [p0 + 1.e-3*np.random.randn(N_dim)*p0 for i in xrange(N_walker)]

    # initialize the EnsembleSampler
    sampler = emcee.EnsembleSampler(Nwalker, Ndim,
                                    fun_prob, args=data,
                                    threads=N_threads)

    # MCMC burn-in run
    pos, prob, state = sampler.run_mcmc(p0, N_burnin)
    sampler.reset()
    fig_res = plt.figure()
    plt.plot(sampler.chain[:, :, 0].T, '-', color='k')
    fig_res.savefig(output_parameter_check)

    # MCMC formal run
    pos, prob, state = sampler.run_mcmc(pos, N_use)
    fig_cnr = triangle.corner(sampler.flatchain,
                              labels=param_labels,
                              truths=np.median(sampler.flatchain[:], axis=0))
    tmp.savefig(output_corner)

    # save chain data
    np.savetxt(output_chaindata, sampler.flatchain[:], delimiter=',')
