# -*- coding: utf-8 -*-
"""

Author
------
Bo Zhang

Email
-----
bozhang@nao.cas.cn

Created on
----------
- Sun Jan  4 20:08:23 2016

Modifications
-------------
-

Aims
----
- implement handy functions used for running emcee

"""


import numpy as np
import emcee
import scipy.optimize as opt
import matplotlib.pyplot as plt
import corner


def emcee_general_run(lnprob,
                      data,
                      param_labels,
                      output_prefix='./mcmc_run',
                      theta0='scipyopt',
                      N=(50, 3, 1, 5000, 1000)):
    """ use this to perform a general MCMC/emcee run

    Parameters:
    -----------

    ln_prob: function
        The function returns the likelihood

    data: array or list of arrays
        The list contains data. This should fit the fun_model.

    param_labels: list of strings
        The list of the parameter names.

    output_prefix: string
        The prefix of output files, used to determine path of outputs.

    theta0: {tuple/list/array or 'scipyopt',}
        Initial guess for theta0.
        The default value is 'scipyopt', which uses scipy.opt function to
        generate a best fit theta0.

    N: tuple (N_walker, N_dim, N_threads, N_burnin, N_use)
        The set of integer numbers used in the MCMC run.
        The default value is (50, 3, 1, 500, 100), whenever fits the data
        or not.
    """
    # unpack parameters
    N_walker, N_dim, N_threads, N_burnin, N_use = N
    output_parameter_check_burnin = \
        u'%s_burnin_%d_use_%d_pcheck_burnin.png' % \
        (output_prefix, N_burnin, N_use)
    output_parameter_check_use = \
        u'%s_burnin_%d_use_%d_pcheck_use.png' % \
        (output_prefix, N_burnin, N_use)
    output_corner = \
        u'%s_burnin_%d_use_%d_corner.png' % \
        (output_prefix, N_burnin, N_use)
    output_chaindata = \
        u'%s_burnin_%d_use_%d_chaindata.csv' % \
        (output_prefix, N_burnin, N_use)

    # initial guess
    if theta0 is 'scipyopt':
        nll = lambda *args: - lnprob(*args)
        p0 = [1., 1., 1., 1.]
        result = opt.minimize(nll, p0, args=data)
        p0 = result['x']
    else:
        p0 = theta0
    p0 = [p0 + 1.e-3*np.random.randn(N_dim) for i in xrange(N_walker)]

    # initialize the EnsembleSampler
    sampler = emcee.EnsembleSampler(N_walker, N_dim,
                                    lnprob, args=data,
                                    threads=N_threads)

    # MCMC burn-in run
    pos, prob, state = sampler.run_mcmc(p0, N_burnin)
    fig_res = plt.figure()
    ax = fig_res.add_subplot(111)
    ax.plot(sampler.chain[:, :, 0].T, '-', color='k')
    fig_res.savefig(output_parameter_check_burnin)
    del fig_res

    # MCMC formal run
    sampler.reset()
    pos, prob, state = sampler.run_mcmc(pos, N_use)
    fig_res = plt.figure()
    ax = fig_res.add_subplot(111)
    ax.plot(sampler.chain[:, :, 0].T, '-', color='k')
    fig_res.savefig(output_parameter_check_use)

    # draw corner
    fig_cnr = corner.corner(sampler.flatchain,
                            labels=param_labels,
                            truths=np.median(sampler.flatchain[:], axis=0))
    fig_cnr.savefig(output_corner)

    # save chain data
    np.savetxt(output_chaindata, sampler.flatchain[:], delimiter=',')

    # return [16., 50., 84.] percentiles of estimated parameters
    return np.percentile(sampler.flatchain[:], [16., 50., 84.], axis=0)


if __name__ == "__main__":
    print "@Cham: test function not implemented ..."
    print "@Cham: but emcee_general_run was tested to be OK ..."