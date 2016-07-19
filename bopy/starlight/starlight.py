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
- Sat Jul 16 22:30:00 2016

Modifications
-------------
- Sat Jul 16 22:30:00 2016    framework
- Sun Jul 17 23:00:00 2016    StarlightGrid
- Mon Jul 18 09:27:00 2016

Aims
----
- to generate config files of STARLIGHT

"""

from __future__ import print_function

import os
import numpy as np
from astropy.table import Table, Column
from astropy.io import fits

__extra_comments__ = '''
# This grid file is written using bopy.starlight
# BOPY, Bo Zhang's python package, https://github.com/hypergravity/bopy
'''


# ################
# Starlight Grid #
# ################


class StarlightGrid(object):
    """ StarlightGrid class is to represent the Grid file object for STARLIGHT
    """
    # specified order of meta data
    meta_order = ['num_of_fits_to_run',
                  'base_dir',
                  'obs_dir',
                  'mask_dir',
                  'out_dir',
                  'phone_number',
                  'llow_SN',
                  'lupp_SN',
                  'Olsyn_ini',
                  'Olsyn_fin',
                  'Odlsyn',
                  'fscale_chi2',
                  'fit_fxk',
                  'IsErrSpecAvailable',
                  'IsFlagSpecAvailable']
    # default values for StarlightGrid instance
    meta = dict(num_of_fits_to_run=2,
                base_dir='/pool/projects/starlight/STARLIGHTv04/BasesDir/',
                obs_dir='/pool/projects/starlight/STARLIGHTv04/',
                mask_dir='/pool/projects/starlight/STARLIGHTv04/',
                out_dir='/pool/projects/starlight/STARLIGHTv04/',
                phone_number='-2007200',
                llow_SN=4730.0,
                lupp_SN=4780.0,
                Olsyn_ini=3400.0,
                Olsyn_fin=8900.0,
                Odlsyn=1.0,
                fit_fxk='FIT',
                fscale_chi2=1.0,
                IsErrSpecAvailable=1,
                IsFlagSpecAvailable=1)
    # specified order of arq
    arq_order = ['arq_obs',
                 'arq_config',
                 'arq_base',
                 'arq_masks',
                 'red_law_option',
                 'v0_start',
                 'vd_start',
                 'arq_out']
    # default arq data
    arq_obs = []
    arq_config = []
    arq_base = []
    arq_masks = []
    red_law_option = []
    v0_start = []
    vd_start = []
    arq_out = []
    # extra comments
    extra = __extra_comments__

    def __init__(self, **kwargs):
        """ initialize instance using arq data """
        for key in kwargs.keys():
            self.__setattr__(key, kwargs[key])

    def set_meta(self, **kwargs):
        """ set meta data """
        for key in kwargs.keys():
            self.meta[key] = kwargs[key]

    def is_arq_validate(self):
        """ check if the arq data length validate """
        try:
            n_obs = len(self.arq_obs)
            assert n_obs == len(self.arq_config) \
                or np.isscalar(self.arq_config)
            assert n_obs == len(self.arq_base) \
                or np.isscalar(self.arq_base)
            assert n_obs == len(self.arq_masks) \
                or np.isscalar(self.arq_masks)
            assert n_obs == len(self.red_law_option) \
                or np.isscalar(self.red_law_option)
            assert n_obs == len(self.v0_start) \
                or np.isscalar(self.v0_start)
            assert n_obs == len(self.vd_start) \
                or np.isscalar(self.vd_start)
            assert n_obs == len(self.arq_out)
            return True
        except AssertionError:
            return False

    def pprint_meta(self):
        """ print meta data """
        for key in self.meta_order:
            print("%s: %s" % (key, self.meta[key]))

    def _print_arq_(self, arq_key):
        """ print single arq field data """
        assert arq_key in self.arq_order

        arq_val = self.__getattribute__(arq_key)

        if np.isscalar(arq_val) or len(arq_val) <= 3:
            print(arq_val)
        else:
            print('[%s, %s, ... %s]' % (arq_val[0], arq_val[1], arq_val[-1]))

    def pprint_arq(self):
        """ print all arq data"""
        for key in self.arq_order:
            print("%s:" % key),
            self._print_arq_(key)

    def pprint(self):
        """ print a summary of the instance """
        print("")
        print('StarlightGrid summary:')
        print('############### [meta] ###############')
        self.pprint_meta()
        print('############### [arq]  ###############')
        self.pprint_arq()
        print(self.extra)
        print('######################################')

    def _meta_to_string(self, val_width=50):
        """ convert meta data to string list """
        fmt_str = "%%%ds" % (-val_width)
        return [(fmt_str + "[%s]\n") % (self.meta[key], key)
                for key in self.meta_order]

    def _arq_to_string(self, sep='   '):
        """ convert arq data to string list """
        #  1. arq data to list
        n_obs = len(self.arq_obs)
        for key in self.arq_order:
            val = self.__getattribute__(key)
            if np.isscalar(val):
                self.__setattr__(key, [val for _ in range(n_obs)])
            else:
                assert len(val) == n_obs

        #  2. to string
        str_list = []
        for i in range(n_obs):
            arq_data_i = ["%s" % self.__getattribute__(key)[i]
                          for key in self.arq_order]
            str_list.append(sep.join(arq_data_i) + "\n")
        return str_list

    def write(self, filepath, meta_val_width=50, sep='   '):
        f = open(filepath, "w+")
        f.writelines(self._meta_to_string(meta_val_width))
        f.write('\n')
        f.writelines(self._arq_to_string(sep=sep))
        f.write('\n')
        f.write(self.extra)
        f.close()


def _test_starlight_grid():
    sg = StarlightGrid(arq_obs=['0414.51901.393.cxt',
                                '0784.52327.478.cxt'],
                       arq_config='StCv04.C99.config',
                       arq_base='Base.BC03.N',
                       arq_masks=['Mask.0414.51901.393.cxt.sc1.CRAP.gm.BN',
                                  'Mask.0784.52327.478.cxt.sc2.CRAP.gm.BN'],
                       red_law_option='CCM',
                       v0_start=0.0,
                       vd_start=150.0,
                       arq_out=['0414.51901.393.cxt.sc4.C99.im.CCM.BN',
                                '0784.52327.478.cxt.sc4.C99.im.CCM.BN']
                       )
    sg.pprint_meta()
    sg.pprint_arq()
    sg.pprint()
    for s in sg._meta_to_string():
        print(s)
    for s in sg._arq_to_string(',:'):
        print(s)
    sg.write('/pool/projects/starlight/STARLIGHTv04/grid_example1.in_')


# ################
# Starlight Base #
# ################


class StarlightBase(object):
    """ StarlightBase class is to represent the Base file object for STARLIGHT
    """
    # specified order of meta data
    meta_order = ['n_base']
    # default values for StarlightBase instance
    meta = dict(n_base=45)
    # specified order of arq data
    arq_order = ['spec_file',
                 'age',
                 'z',
                 'code',
                 'mstar',
                 'yav',
                 'afe']
    # default values for bases
    spec_file = []
    age = []
    z = []
    code = []
    mstar = []
    yav = []
    afe = []
    # extra comments
    extra = __extra_comments__

    def __init__(self, **kwargs):
        """ initialize instance using arq data """
        for key in kwargs.keys():
            self.__setattr__(key, kwargs[key])
        # count spec_file
        self.meta['n_base'] = len(self.spec_file)

    def _meta_to_string(self, val_width=50):
        """ convert meta data to string list """
        fmt_str = "%%%ds" % (-val_width)
        return [(fmt_str + "[%s]\n") % (self.meta[key], key)
                for key in self.meta_order]

    def _arq_to_string(self, sep='   '):
        """ convert arq data to string list """
        #  1. arq data to list
        n_obs = len(self.spec_file)
        for key in self.arq_order:
            val = self.__getattribute__(key)
            if np.isscalar(val):
                self.__setattr__(key, [val for _ in range(n_obs)])
            else:
                assert len(val) == n_obs

        # 2. to string
        str_list = []
        for i in range(n_obs):
            arq_data_i = ["%s" % self.__getattribute__(key)[i]
                          for key in self.arq_order]
            str_list.append(sep.join(arq_data_i) + "\n")
        return str_list

    def write(self, filepath, meta_val_width=50, sep='   '):
        f = open(filepath, "w+")
        f.writelines(self._meta_to_string(meta_val_width))
        # f.write('\n')
        f.writelines(self._arq_to_string('   '))
        f.write('\n')
        f.write(self.extra)
        f.close()


def _test_starlight_base():
    sg = StarlightBase(spec_file=['bc2003_hr_m42_chab_ssp_020.spec',
                                  'bc2003_hr_m42_chab_ssp_045.spec'],
                       age=[0.00100e9,
                            0.00316e9],
                       z=[0.00400,
                          0.00400],
                       code=['age020_m42',
                             'age020_m42'],
                       mstar=[1.0000,
                              0.9999],
                       yav=[0,
                            0],
                       afe=[0.0,
                            0.0]
                       )
    sg.write('/pool/projects/starlight/STARLIGHTv04/Base.BC03.N_')


_config_StCv04C99_ = dict(
    # Normalization lambdas
    l_norm=4020.0,
    llow_norm=4010.0,
    lupp_norm=4060.0,
    # Parameter Limits
    AV_low=-1.0,
    AV_upp=4.0,
    YAV_low=-0.0001,
    YAV_upp=0.0001,
    fn_low=0.7,
    fn_upp=1.3,
    v0_low=-500.0,
    v0_upp=500.0,
    vd_low=0.0,
    vd_upp=500.0,
    # Clipping options & Weight-Control-Filter
    clip_method_option='NSIGMA', # NOCLIP/NSIGMA/RELRES/ABSRES
    sig_clip_threshold=3.0,
    wei_nsig_threshold=2.0,
    # Miscellaneous
    dl_cushion=50.0,
    f_cut=0.001,
    N_int_Gauss=31,
    i_verbose=1,
    i_verbose_anneal=1,
    Is1stLineHeader=0,
    i_FastBC03_FLAG=1,
    i_FitPowerLaw=0,
    alpha_PowerLaw=-0.5,
    i_SkipExistingOutFiles=0,
    # Markov Chains technical parameters
    N_chains=7,
    xinit_max=0.50,
    i_UpdateEps=0,
    i_UpdateAlpha=2,
    Falpha=2.0,
    i_MoveOneParOnly=1,
    i_UpdateAVYAVStepSeparately=1,
    i_HelpParWithMaxR=1,
    prob_jRmax=0.2,
    i_HelpPopVectorMove2Average=1,
    prob_HelpPopVectorMove2Average=0.4,
    i_HelpAVYAVMove2Average=1,
    prob_HelpAVYAVMove2Average=0.4,
    NRC_AV_Default=10,
    # First Fits (FF) technical parameters
    Temp_ini_FF=1.0e2,
    Temp_fin_FF=1.0,
    fini_eps_FF=1.0e1,
    ffin_eps_FF=1.0e2,
    R_ini_FF=1.3,
    R_fin_FF=1.2,
    IsGRTestHard_FF=0,
    N_loops_FF=10,
    i_RestartChains_FF=1,
    fN_sim_FF=1.0e1,
    fNmax_steps_FF=1.0e4,
    eff_IDEAL_FF=0.23,
    # GR R-threshold & Method for Burn-In loop
    R_Burn_in=1.2,
    IsGRTestHard_BurnIn=1,
    # EX0s technical parameters
    EXOs_PopVector_option='MIN', # MIN/AVE
    EXOs_method_option='CUMUL', # CUMUL/SMALL
    EXOs_Threshold=0.02,
    Temp_ini_EX0s=1.0,
    Temp_fin_EX0s=1.0e-3,
    fini_eps_EX0s=1.0e2,
    ffin_eps_EX0s=1.0e3,
    R_ini_EX0s=1.2,
    R_fin_EX0s=1.0,
    IsGRTestHard_EX0s=1,
    N_loops_EX0s=10,
    i_RestartChains_EX0s=1,
    fN_sim_EX0s=1.0e2,
    fNmax_steps_EX0s=1.0e3,
    eff_IDEAL_EX0s=0.50,
    IsScaleNstepsInEX0sFits=1,
    IsNoKin4LargeBaseInEX0sFits=0,
    frac_NoKin4LargeBaseInEX0sFits=0.0,
    fEX0_MinBaseSize=0.1)
_config_StCv04C11_ = dict(
    # Normalization lambdas
    l_norm=4020.0,
    llow_norm=4010.0,
    lupp_norm=4060.0,
    # Parameter Limits
    AV_low=-1.0,
    AV_upp=4.0,
    YAV_low=-0.0001,
    YAV_upp=0.0001,
    fn_low=0.7,
    fn_upp=1.3,
    v0_low=-500.0,
    v0_upp=500.0,
    vd_low=0.0,
    vd_upp=500.0,
    # Clipping options & Weight-Control-Filter
    clip_method_option='NSIGMA', # NOCLIP/NSIGMA/RELRES/ABSRES
    sig_clip_threshold=3.0,
    wei_nsig_threshold=2.0,
    # Miscellaneous
    dl_cushion=50.0,
    f_cut=0.001,
    N_int_Gauss=31,
    i_verbose=1,
    i_verbose_anneal=0,
    Is1stLineHeader=0,
    i_FastBC03_FLAG=1,
    i_FitPowerLaw=0,
    alpha_PowerLaw=-0.5,
    i_SkipExistingOutFiles=0,
    # Markov Chains technical parameters
    N_chains=7,
    xinit_max=0.50,
    i_UpdateEps=0,
    i_UpdateAlpha=2,
    Falpha=2.0,
    i_MoveOneParOnly=1,
    i_UpdateAVYAVStepSeparately=1,
    i_HelpParWithMaxR=1,
    prob_jRmax=0.2,
    i_HelpPopVectorMove2Average=1,
    prob_HelpPopVectorMove2Average=0.4,
    i_HelpAVYAVMove2Average=1,
    prob_HelpAVYAVMove2Average=0.4,
    NRC_AV_Default=10,
    # First Fits (FF) technical parameters
    Temp_ini_FF=1.0e2,
    Temp_fin_FF=1.0,
    fini_eps_FF=1.0e1,
    ffin_eps_FF=1.0e2,
    R_ini_FF=1.3,
    R_fin_FF=1.3,
    IsGRTestHard_FF=0,
    N_loops_FF=3,
    i_RestartChains_FF=1,
    fN_sim_FF=1.0e1,
    fNmax_steps_FF=1.0e4,
    eff_IDEAL_FF=0.23,
    # GR R-threshold & Method for Burn-In loop
    R_Burn_in=1.2,
    IsGRTestHard_BurnIn=0,
    # EX0s technical parameters
    EXOs_PopVector_option='MIN',  # MIN/AVE
    EXOs_method_option='CUMUL',  # CUMUL/SMALL
    EXOs_Threshold=0.02,
    Temp_ini_EX0s=1.0,
    Temp_fin_EX0s=1.0e-3,
    fini_eps_EX0s=1.0e2,
    ffin_eps_EX0s=1.0e3,
    R_ini_EX0s=1.2,
    R_fin_EX0s=1.0,
    IsGRTestHard_EX0s=1,
    N_loops_EX0s=5,
    i_RestartChains_EX0s=1,
    fN_sim_EX0s=1.0e2,
    fNmax_steps_EX0s=1.0e3,
    eff_IDEAL_EX0s=0.50,
    IsScaleNstepsInEX0sFits=1,
    IsNoKin4LargeBaseInEX0sFits=1,
    frac_NoKin4LargeBaseInEX0sFits=0.0,
    fEX0_MinBaseSize=0.1)

_config_comments_ = (
    ('# Configuration parameters for StarlightChains_v04.for - Cid@Lagoa - 18/Feb/2007 #\n'
     '#\n#\n# Normalization lambdas\n#\n'),
    '#\n#\n# Parameter Limits\n#\n',
    '#\n#\n# Clipping options & Weight-Control-Filter\n#\n',
    '#\n#\n# Miscellaneous\n#\n',
    '#\n#\n# Markov Chains technical parameters\n#\n',
    '#\n#\n# First Fits (FF) technical parameters\n#\n',
    '#\n#\n# GR R-threshold & Method for Burn-In loop\n#\n',
    '#\n#\n# EX0s technical parameters\n#\n',
    (
        '\n\nCid@Lagoa - 18/February/2007'
        '\n\nOBS: A SLOW config!'
        '\n\n\n\n'
        'Technical parameters you may want to play with to obtain FAST, MEDIUM'
        '\n& SLOW fits:\n\n'
        '--------------------------------\n'
        '|  FAST  |    MEDIUM  |  LONG  |\n'
        '--------------------------------\n'
        '|  5     |    7       |  12    | [N_chains]\n'
        '|  3     |    5       |  10    | [N_loops_FF & *_EX0s]\n'
        '|  1.3   |    1.2     |  1.1   | [R_ini_FF & R_fin_FF & *_EX0s]\n'
        '|  0     |   0 or 1   |  1     | [IsGRTestHard_FF & *_BurnIn & *_EX0s]\n'
        '--------------------------------\n')
    )
_config_comments_insert_index_ = (
    0, 4, 15, 19, 30, 45, 58, 61, 65)


# ##################
# Starlight Config #
# ##################


class StarlightConfig(object):
    """ StarlightConfig class is to represent the Config file object for STARLIGHT
    """
    # specified order of meta data
    meta_order = [
        # Normalization lambdas
        'l_norm',
        'llow_norm',
        'lupp_norm',
        # Parameter Limits
        'AV_low',
        'AV_upp',
        'YAV_low',
        'YAV_upp',
        'fn_low',
        'fn_upp',
        'v0_low',
        'v0_upp',
        'vd_low',
        'vd_upp',
        # Clipping options & Weight-Control-Filter
        'clip_method_option',
        'sig_clip_threshold',
        'wei_nsig_threshold',
        # Miscellaneous
        'dl_cushion',
        'f_cut',
        'N_int_Gauss',
        'i_verbose',
        'i_verbose_anneal',
        'Is1stLineHeader',
        'i_FastBC03_FLAG',
        'i_FitPowerLaw',
        'alpha_PowerLaw',
        'i_SkipExistingOutFiles',
        # Markov Chains technical parameters
        'N_chains',
        'xinit_max',
        'i_UpdateEps',
        'i_UpdateAlpha',
        'Falpha',
        'i_MoveOneParOnly',
        'i_UpdateAVYAVStepSeparately',
        'i_HelpParWithMaxR',
        'prob_jRmax',
        'i_HelpPopVectorMove2Average',
        'prob_HelpPopVectorMove2Average',
        'i_HelpAVYAVMove2Average',
        'prob_HelpAVYAVMove2Average',
        'NRC_AV_Default',
        # First Fits (FF) technical parameters
        'Temp_ini_FF',
        'Temp_fin_FF',
        'fini_eps_FF',
        'ffin_eps_FF',
        'R_ini_FF',
        'R_fin_FF',
        'IsGRTestHard_FF',
        'N_loops_FF',
        'i_RestartChains_FF',
        'fN_sim_FF',
        'fNmax_steps_FF',
        'eff_IDEAL_FF',
        # GR R-threshold & Method for Burn-In loop
        'R_Burn_in',
        'IsGRTestHard_BurnIn',
        # EX0s technical parameters
        'EXOs_PopVector_option',
        'EXOs_method_option',
        'EXOs_Threshold',
        'Temp_ini_EX0s',
        'Temp_fin_EX0s',
        'fini_eps_EX0s',
        'ffin_eps_EX0s',
        'R_ini_EX0s',
        'R_fin_EX0s',
        'IsGRTestHard_EX0s',
        'N_loops_EX0s',
        'i_RestartChains_EX0s',
        'fN_sim_EX0s',
        'fNmax_steps_EX0s',
        'eff_IDEAL_EX0s',
        'IsScaleNstepsInEX0sFits',
        'IsNoKin4LargeBaseInEX0sFits',
        'frac_NoKin4LargeBaseInEX0sFits',
        'fEX0_MinBaseSize']
    # default values for StarlightConfig instance (StCv04.C99.config)
    meta = _config_StCv04C99_
    # extra comments
    extra = __extra_comments__
    # necessary comments
    config_comments = _config_comments_
    config_comments_insert_index = _config_comments_insert_index_

    def __init__(self, **kwargs):
        """ initialize StarlightConfig instance with StCv04.C99.config """
        self.set_meta(**kwargs)

    def set_meta(self, **kwargs):
        """ set meta data """
        for key in kwargs.keys():
            self.meta[key] = kwargs[key]

    def set_quick(self, template='StCv04.C99.config'):
        if template == 'StCv04.C99.config':
            self.set_meta(**_config_StCv04C99_)
        elif template == 'StCv04.C11.config':
            self.set_meta(**_config_StCv04C11_)
        else:
            raise(ValueError('@Cham: your option should be one of {''StCv04.C99.config'', ''StCv04.C11.config''}'))

    def _meta_to_string(self, val_width=50):
        """ convert meta data to string list """
        fmt_str = "%%%ds" % (-val_width)
        return [(fmt_str + "[%s]\n") % (self.meta[key], key)
                for key in self.meta_order]

    def write(self, filepath, meta_val_width=50):
        f = open(filepath, "w+")
        str_list = self._meta_to_string(meta_val_width)

        # insert necessary comments
        for i in range(len(self.config_comments_insert_index)-1):
            str_list.insert(
                self.config_comments_insert_index[i], self.config_comments[i])
        str_list.append(self.config_comments[-1])

        # for str in str_list:
        #     print(str)

        f.writelines(str_list)
        f.write('\n')
        f.write(self.extra)
        f.close()


def _test_starlight_config():
    sc = StarlightConfig()
    sc.set_quick('StCv04.C99.config')
    print(sc.meta)
    sc.set_quick('StCv04.C11.config')
    print(sc.meta)
    sc.set_meta(Is1stLineHeader=1)
    sc.write('/pool/projects/starlight/STARLIGHTv04/StCv04.C99.config_')


# ################
# Starlight Mask #
# ################
sm_tplt_sdss_gm = Table(np.array([
    [3710.000000, 3744.000000, 0.000000,
     '[OII]          3726+3729 emission line', 0],
    [3858.000000, 3880.000000, 0.000000,
     '[NeIII]        3869 emission line', 0],
    [3960.000000, 3980.000000, 0.000000,
     'Hepsilon       3970 emission line (over CaII H)', 0],
    [4092.000000, 4112.000000, 0.000000,
     'Hdelta         4102 emission line', 0],
    [4330.000000, 4350.000000, 0.000000,
     'Hgamma         4340 emission line', 0],
    [4848.000000, 4874.000000, 0.000000,
     'Hbeta          4861 emission line', 0],
    [4940.000000, 5028.000000, 0.000000,
     '[OIII]         4959 & 5007 emission lines', 0],
    [5866.000000, 5916.000000, 0.000000,
     'HeI & NaD      5876 emission line & ~ 5890 ISM absorption', 0],
    [6280.000000, 6320.000000, 0.000000,
     '[OI]           6300 emission line', 0],
    [6528.000000, 6608.000000, 0.000000,
     'Halpha & [NII] 6563, 6548 & 6583 emission lines', 0],
    [6696.000000, 6752.000000, 0.000000,
     '[SII]          6717 & 6731 emission lines', 0]]),
    names=['wave1', 'wave2', 'mask_value', 'comment', 'sky'],
    dtype=['f8', 'f8', 'f8', 'S70', '?'])
sm_tplt_bc03_long = Table(np.array([
    [3718.000000, 3735.000000, 0.000000, '*[OII]     3726.0', 0],
    [3720.000000, 3737.000000, 0.000000, '*[OII]     3728.8', 0],
    [3796.000000, 3802.000000, 0.000000, 'H10        3797.9', 0],
    [3833.000000, 3841.000000, 0.000000, 'H90        3835.4', 0],
    [3862.000000, 3871.000000, 0.000000, '*[NeIII]   3869.1', 0],
    [3882.000000, 3891.000000, 0.000000, 'H8HeI      3889.0', 0],
    [3961.000000, 3973.000000, 0.000000, '[NeIII]    3967.8', 0],
    [3961.000000, 3974.000000, 0.000000, 'Hepsilon   3970.1', 0],
    [4024.000000, 4030.000000, 0.000000, 'HeI        4026.2', 0],
    [4067.000000, 4070.000000, 0.000000, '[SII]      4068.6', 0],
    [4095.000000, 4107.000000, 0.000000, '*Hdelta    4101.7', 0],
    [4332.000000, 4346.000000, 0.000000, '*Hgamma    4340.5', 0],
    [4469.000000, 4475.000000, 0.000000, 'HeI        4471.5', 0],
    [4854.000000, 4868.000000, 0.000000, '*Hbeta     4861.3', 0],
    [4947.000000, 4968.000000, 0.000000, '*[OIII]    4958.9', 0],
    [4998.000000, 5019.000000, 0.000000, '*[OIII]    5006.8', 0],
    [5009.000000, 5022.000000, 0.000000, 'HeI        5015.7', 0],
    [5156.000000, 5164.000000, 0.000000, '[FeVII]    5158.4', 0],
    [5193.000000, 5202.000000, 0.000000, '[NI]       5199.1', 0],
    [5267.000000, 5274.000000, 0.000000, '[FeIII]    5270.4', 0],
    [5868.000000, 5883.000000, 0.000000, 'HeI        5875.6', 0],
    [6296.000000, 6303.000000, 0.000000, '*[OI]      6300.3', 0],
    [6303.000000, 6315.000000, 0.000000, '[SIII]     6312.1', 0],
    [6355.000000, 6372.000000, 0.000000, '[OI]       6363.8', 0],
    [6538.000000, 6591.000000, 0.000000, '*[NII]     6548.0', 0],
    [6547.000000, 6576.000000, 0.000000, '*Halpha    6562.8', 0],
    [6575.000000, 6591.000000, 0.000000, '*[NII]     6583.5', 0],
    [6667.000000, 6681.000000, 0.000000, 'HeI        6678.2', 0],
    [6710.000000, 6743.000000, 0.000000, '*[SII]     6716.4', 0],
    [6718.000000, 6743.000000, 0.000000, '*[SII]     6730.8', 0],
    [7061.000000, 7072.000000, 0.000000, 'HeI        7065.2', 0],
    [7125.000000, 7142.000000, 0.000000, '[ArIII]    7135.8', 0],
    [7315.000000, 7329.000000, 0.000000, '[OII]      7319.5', 0],
    [7320.000000, 7335.000000, 0.000000, '[OII]      7330.2', 0],
    [3718.000000, 3735.000000, 0.000000, '*[OII]     3726.0', 0],
    [3720.000000, 3737.000000, 0.000000, '*[OII]     3728.8', 0],
    [3796.000000, 3802.000000, 0.000000, 'H10        3797.9', 0],
    [3833.000000, 3841.000000, 0.000000, 'H90        3835.4', 0],
    [3862.000000, 3871.000000, 0.000000, '*[NeIII]   3869.1', 0],
    [3882.000000, 3891.000000, 0.000000, 'H8HeI      3889.0', 0],
    [3961.000000, 3973.000000, 0.000000, '[NeIII]    3967.8', 0],
    [3961.000000, 3974.000000, 0.000000, 'Hepsilon   3970.1', 0],
    [4024.000000, 4030.000000, 0.000000, 'HeI        4026.2', 0],
    [4067.000000, 4070.000000, 0.000000, '[SII]      4068.6', 0],
    [4095.000000, 4107.000000, 0.000000, '*Hdelta    4101.7', 0],
    [4332.000000, 4346.000000, 0.000000, '*Hgamma    4340.5', 0],
    [4469.000000, 4475.000000, 0.000000, 'HeI        4471.5', 0],
    [4854.000000, 4868.000000, 0.000000, '*Hbeta     4861.3', 0],
    [4947.000000, 4968.000000, 0.000000, '*[OIII]    4958.9', 0],
    [4998.000000, 5019.000000, 0.000000, '*[OIII]    5006.8', 0],
    [5009.000000, 5022.000000, 0.000000, 'HeI        5015.7', 0],
    [5156.000000, 5164.000000, 0.000000, '[FeVII]    5158.4', 0],
    [5193.000000, 5202.000000, 0.000000, '[NI]       5199.1', 0],
    [5267.000000, 5274.000000, 0.000000, '[FeIII]    5270.4', 0],
    [5868.000000, 5883.000000, 0.000000, 'HeI        5875.6', 0],
    [6296.000000, 6303.000000, 0.000000, '*[OI]      6300.3', 0],
    [6303.000000, 6315.000000, 0.000000, '[SIII]     6312.1', 0],
    [6355.000000, 6372.000000, 0.000000, '[OI]       6363.8', 0],
    [6538.000000, 6591.000000, 0.000000, '*[NII]     6548.0', 0],
    [6547.000000, 6576.000000, 0.000000, '*Halpha    6562.8', 0],
    [6575.000000, 6591.000000, 0.000000, '*[NII]     6583.5', 0],
    [6667.000000, 6681.000000, 0.000000, 'HeI        6678.2', 0],
    [6710.000000, 6743.000000, 0.000000, '*[SII]     6716.4', 0],
    [6718.000000, 6743.000000, 0.000000, '*[SII]     6730.8', 0],
    [7061.000000, 7072.000000, 0.000000, 'HeI        7065.2', 0],
    [7125.000000, 7142.000000, 0.000000, '[ArIII]    7135.8', 0],
    [7315.000000, 7329.000000, 0.000000, '[OII]      7319.5', 0],
    [7320.000000, 7335.000000, 0.000000, '[OII]      7330.2', 0],
    [6845.000000, 6945.000000, 0.000000,
     'BC03-bug   shifted -5 Angs wrt to what they say.', 0],
    [7165.000000, 7210.000000, 0.000000,
     'BC03-bug?  looks like a bug, but not mentioned in BC03.', 0],
    [7550.000000, 7725.000000, 0.000000,
     'BC03-bug   equal to what they say in their paper.', 0],
    [5880.000000, 5906.000000, 0.000000, 'NaD 5890 & 5896 (gm!)', 0],
    [9059.000000, 9079.000000, 0.000000, '[SIII]     9069 (gm!)', 0]]),
    names=['wave1', 'wave2', 'mask_value', 'comment', 'sky'],
    dtype=['f8', 'f8', 'f8', 'S70', '?'])
sm_tplt_bc03_short = Table(np.array([
    [3721.000000, 3728.000000, 0.000000, '*[OII]     3726.0', 0],
    [3721.000000, 3730.000000, 0.000000, '*[OII]     3728.8', 0],
    [3746.000000, 3749.000000, 0.000000, 'H12        3750.2', 0],
    [3768.000000, 3773.000000, 0.000000, 'H11        3770.6', 0],
    [3887.000000, 3892.000000, 0.000000, 'H8HeI      3889.0', 0],
    [4023.000000, 4025.000000, 0.000000, 'HeI        4026.2', 0],
    [4065.000000, 4076.000000, 0.000000, '[SII]      4068.6', 0],
    [4654.000000, 4658.000000, 0.000000, '[FeIII]    4658.0', 0],
    [4915.000000, 4921.000000, 0.000000, 'HeI        4921.9', 0],
    [4951.000000, 4962.000000, 0.000000, '*[OIII]    4958.9', 0],
    [4999.000000, 5010.000000, 0.000000, '*[OIII]    5006.8', 0],
    [5262.000000, 5274.000000, 0.000000, '[FeIII]    5270.4', 0],
    [6082.000000, 6085.000000, 0.000000, '[FeVII]    6086.3', 0],
    [6305.000000, 6320.000000, 0.000000, '[SIII]     6312.1', 0],
    [6539.000000, 6561.000000, 0.000000, '*[NII]     6548.0', 0],
    [6574.000000, 6596.000000, 0.000000, '*[NII]     6583.5', 0],
    [6673.000000, 6677.000000, 0.000000, 'HeI        6678.2', 0],
    [6721.000000, 6744.000000, 0.000000, '*[SII]     6730.8', 0],
    [6845.000000, 6945.000000, 0.000000,
     'BC03-bug   shifted -5 Angs wrt to what they say.', 0],
    [7165.000000, 7210.000000, 0.000000,
     'BC03-bug?  looks like a bug, but not mentioned in BC03.', 0],
    [7550.000000, 7725.000000, 0.000000,
     'BC03-bug   equal to what they say in their paper.', 0],
    [5880.000000, 5906.000000, 0.000000, 'NaD 5890 & 5896 (gm!)', 0],
    [9059.000000, 9079.000000, 0.000000, '[SIII]     9069 (gm!)', 0]]),
    names=['wave1', 'wave2', 'mask_value', 'comment', 'sky'],
    dtype=['f8', 'f8', 'f8', 'S70', '?'])


class StarlightMask(Table):
    """ StarlightMask class is to represent the Mask file object for STARLIGHT
    """

    def __init__(self):
        super(self.__class__, self).__init__(
            data=[[], [], [], [], []],
            names=['wave1', 'wave2', 'mask_value', 'comment', 'sky'],
            dtype=['f8', 'f8', 'f8', 'S70', '?'],
            masked=True)

    def quick_set(self, template='sdss_gm'):
        """ a quick set to Starlight Mask template """
        # options
        map_dict = {'sdss_gm': sm_tplt_sdss_gm,
                    'bc03_short': sm_tplt_bc03_short,
                    'bc03_long': sm_tplt_bc03_long}
        assert template in map_dict.keys()

        # clean
        self.clean()

        # fill self with template
        for _ in map_dict[template]:
            self.add_row(_)

    def clean(self):
        """ clean all items in a StarlightMask instance """
        while len(self) > 0:
            self.remove_row(0)

    def _to_string(self, sep='  ', z=None):
        """ convert to string list """
        if self['sky'].data.any():
            # need redshift to calculate the wavelength
            assert z is not None

        str_list = []
        for i in range(len(self)):
            if self['sky'][i]:
                str_list.append('%s%s%s%s%s%s%s%s%s\n'
                                % (self['wave1'][i]/(1+z), sep,
                                   self['wave2'][i]/(1+z), sep,
                                   self['mask_value'][i], sep,
                                   self['comment'][i], sep,
                                   self['sky'][i]))
            else:
                str_list.append('%s%s%s%s%s%s%s%s%s\n'
                                % (self['wave1'][i], sep,
                                   self['wave2'][i], sep,
                                   self['mask_value'][i], sep,
                                   self['comment'][i], sep,
                                   self['sky'][i]))
        return str_list

    def write_to(self, filepath, z=None, sep='  '):
        """ write mask in STARLIGHT form """
        f = open(filepath, 'w+')
        f.write('%d\n' % len(self))
        f.writelines(self._to_string(sep=sep, z=z))
        f.write(__extra_comments__)
        f.close()


def _test_starlight_mask():
    sm = StarlightMask()
    sm.quick_set('sdss_gm')
    sm.quick_set('bc03_short')
    sm.add_row([6280, 6288, 2, 'DIB6284', 0])
    sm.add_row([5567, 5587, 0, 'skyline', 1])
    sm.pprint()
    sm.write_to('/pool/projects/starlight/STARLIGHTv04/Masks.EmLines.SDSS.gm_', z=0.01)


# ##################
# Starlight Output #
# ##################


class StarlightOutput(object):
    """ StarlightOutput class is to read/re-construct the STARLIGHT results
    """
    meta = dict(
        arq_obs='',
        arq_base='',
        arq_masks='',
        arq_config='',
        N_base=0,
        N_YAV_components=0,
        i_FitPowerLaw=0,
        alpha_PowerLaw=0,
        red_law_option='',
        q_norm=0,
        l_ini=0.,
        l_fin=0.,
        dl=0.,
        l_norm=0.,
        llow_norm=0.,
        lupp_norm=0.,
        fobs_norm=0.,
        llow_SN=0.,
        lupp_SN=0.,
        SN_in_SN_window=0.,
        SN_in_norm_window=0.,
        SN_err_in_SN_window=0.,
        SN_err_in_norm_window=0.,
        fscale_chi2=0.,
        idum_orig=0,
        NOl_eff=0,
        Nl_eff=0,
        Ntot_cliped=0,
        clip_method='',
        Nglobal_steps=0,
        N_chains=0,
        NEX0s_base=0,
        Clip_Bug=0,
        RC_Crash=0,
        Burn_In_warning_flags=0,
        n_censored_weights=0,
        wei_nsig_threshold=0,
        wei_limit=0,
        idt_all=0,
        wdt_TotTime=0.,
        wdt_UsrTime=0.,
        wdt_SysTime=0.,
        chi2_Nl_eff=0.,
        adev=0.,
        sum_of_x=0.,
        Flux_tot=0.,
        Mini_tot=0.,
        Mcor_tot=0.,
        v0_min=0.,
        vd_min=0.,
        AV_min=0.,
        YAV_min=0.)
    syn_spec = None
    syn_model = None

    def __init__(self, filepath):
        """ initialize instance

        Parameters
        ----------
        filepath: string
            file path of the starlight output

        Returns
        -------
        so: StarlightOutput instance
            StarlightOutput instance

        """
        # assert filepath existence
        try:
            assert os.path.exists(filepath)
        except AssertionError:
            raise(AssertionError('@Cham: file does not exist! %s' % filepath))

        # read header
        f = open(filepath, 'r')
        lines = f.readlines()
        self.meta, state = read_starlight_output_header(lines)
        f.close()
        try:
            assert state
        except AssertionError:
            raise(AssertionError('@Cham: starlight read header FAILED!'))

        # read blocks
        # As pointed out by the STARLIGHT manual, 3 out of 5 blocks are ignored
        # Only the synthetic results (coefficients [1/5] and spectra [5/5])
        # are loaded into data

        # 1. syn model
        N_base = self.meta['N_base']
        self.syn_model = read_starlight_output_syn_model(lines[63:63+N_base])

        # 2. syn spectrum
        syn_spec_start = 63+N_base+5+N_base+2+N_base+11
        Nl_obs = np.int(lines[syn_spec_start-1].split('[')[0].strip())
        # assert that the number of rest lines is equal to Nl_obs
        assert len(lines) - syn_spec_start == Nl_obs
        self.syn_spec = read_starlight_output_syn_spec(lines[syn_spec_start:])

    def pprint(self):
        len_max = np.max([len(key) for key in self.meta.keys()])
        fmt_str = '%%%ds' % (len_max+1)
        for k,v in self.meta.items():
            print((fmt_str+': %s') % (k, v))
        print('')
        print(self.syn_model)
        print('')
        print(self.syn_spec)

    def write_fits(self, filepath, **kwargs):
        prihdr = fits.Header()
        prihdr['AUTHOR'] = 'Bo Zhang'
        for k, v in self.meta.items():
            prihdr[k] = v
        prihdu = fits.PrimaryHDU(header=prihdr)
        hdulist = fits.HDUList([prihdu,
                                fits.table_to_hdu(self.syn_model),
                                fits.table_to_hdu(self.syn_spec)])
        if os.path.exists(filepath):
            print('[StarlightOuput.write_fits()]: filepath exists: %s'
                  % filepath)
        hdulist.writeto(filepath, **kwargs)


def read_starlight_output_syn_spec(lines):
    """ read syn_spec of starlight output """
    Nl_obs = len(lines)
    wave = Column(np.zeros((Nl_obs, ), dtype=np.float), 'wave')
    flux_obs = Column(np.zeros((Nl_obs, ), dtype=np.float), 'flux_obs')
    flux_syn = Column(np.zeros((Nl_obs, ), dtype=np.float), 'flux_syn')
    weight = Column(np.zeros((Nl_obs, ), dtype=np.float), 'weight')
    for i, line in enumerate(lines):
        line_split = line.split()
        wave[i] = np.float(line_split[0])
        flux_obs[i] = np.float(line_split[1])
        flux_syn[i] = np.float(line_split[2])
        weight[i] = np.float(line_split[3])
    return Table([wave, flux_obs, flux_syn, weight])


def read_starlight_output_syn_model(lines):
    """ read syn_model of starlight output """
    N_base=len(lines)
    j = Column(np.zeros((N_base,), dtype=np.int), 'j')
    x_j = Column(np.zeros((N_base,), dtype=np.float), 'x_j')
    Mini_j = Column(np.zeros((N_base,), dtype=np.float), 'Mini_j')
    Mcor_j = Column(np.zeros((N_base,), dtype=np.float), 'Mcor_j')
    age_j = Column(np.zeros((N_base,), dtype=np.float), 'age_j')
    Z_j = Column(np.zeros((N_base,), dtype=np.float), 'Z')
    LM_j = Column(np.zeros((N_base,), dtype=np.float), 'LM_j')
    YAV = Column(np.zeros((N_base,), dtype=np.float), 'YAV')
    Mstars = Column(np.zeros((N_base,), dtype=np.float), 'Mstars')
    component_j = Column(np.zeros((N_base,), dtype=np.string_), 'component_j')
    aFe = Column(np.zeros((N_base,), dtype=np.float), 'aFe')
    SSP_chi2r = Column(np.zeros((N_base,), dtype=np.float), 'SSP_chi2r')
    SSP_adev = Column(np.zeros((N_base,), dtype=np.float), 'SSP_adev')
    SSP_AV = Column(np.zeros((N_base,), dtype=np.float), 'SSP_AV')
    SSP_x = Column(np.zeros((N_base,), dtype=np.float), 'SSP_x')
    for i in range(len(lines)):
        line_split = lines[i].split()
        j[i] = np.int(line_split[0])
        x_j[i] = np.float(line_split[1])
        Mini_j[i] = np.float(line_split[2])
        Mcor_j[i] = np.float(line_split[3])
        age_j[i] = np.float(line_split[4])
        Z_j[i] = np.float(line_split[5])
        LM_j[i] = np.float(line_split[6])
        YAV[i] = np.float(line_split[7])
        Mstars[i] = np.float(line_split[8])
        component_j[i] = line_split[9]
        aFe[i] = np.float(line_split[10])
        SSP_chi2r[i] = np.float(line_split[11])
        SSP_adev[i] = np.float(line_split[12])
        SSP_AV[i] = np.float(line_split[13])
        SSP_x[i] = np.float(line_split[14])
    return Table([j, x_j, Mini_j, Mcor_j, age_j, Z_j, LM_j, YAV, Mstars,
                  component_j, aFe, SSP_chi2r, SSP_adev, SSP_AV, SSP_x])


def read_starlight_output_header(lines):
    """ read header of starlight output
    Although this method is not elegant, it works!
    """
    # initial state False
    state = False
    # initialize meta
    meta = dict()
    try:
        # Some input info
        meta['arq_obs'] = lines[5].split('[')[0].strip()
        meta['arq_base'] = lines[6].split('[')[0].strip()
        meta['arq_masks'] = lines[7].split('[')[0].strip()
        meta['arq_config'] = lines[8].split('[')[0].strip()
        meta['N_base'] = np.int(lines[9].split('[')[0].strip())
        meta['N_YAV_components'] = np.int(lines[10].split('[')[0].strip())
        meta['i_FitPowerLaw'] = np.int(lines[11].split('[')[0].strip())
        meta['alpha_PowerLaw'] = np.float(lines[12].split('[')[0].strip())
        meta['red_law_option'] = lines[13].split('[')[0].strip()
        meta['q_norm'] = np.float(lines[14].split('[')[0].strip())
        # (Re)Sampling Parameters
        meta['l_ini'] = np.float(lines[17].split('[')[0].strip())
        meta['l_fin'] = np.float(lines[18].split('[')[0].strip())
        meta['dl'] = np.float(lines[19].split('[')[0].strip())
        # Normalization info
        meta['l_norm'] = np.float(lines[22].split('[')[0].strip())
        meta['llow_norm'] = np.float(lines[23].split('[')[0].strip())
        meta['lupp_norm'] = np.float(lines[24].split('[')[0].strip())
        meta['fobs_norm'] = np.float(lines[25].split('[')[0].strip())
        # S/N
        meta['llow_SN'] = np.float(lines[28].split('[')[0].strip())
        meta['lupp_SN'] = np.float(lines[29].split('[')[0].strip())
        meta['SN_in_SN_window'] = np.float(lines[30].split('[')[0].strip())
        meta['SN_in_norm_window'] = np.float(lines[31].split('[')[0].strip())
        meta['SN_err_in_SN_window'] = np.float(lines[32].split('[')[0].strip())
        meta['SN_err_in_norm_window'] =\
            np.float(lines[33].split('[')[0].strip())
        meta['fscale_chi2'] = np.float(lines[34].split()[0].strip('['))
        # etc... [ignored ugly data form! --> this makes me happy!]
        meta['NOl_eff'] = np.int(lines[38].split('[')[0].strip())
        meta['Nl_eff'] = np.int(lines[39].split('[')[0].strip())
        lines_40_split = lines[40].split('[')[0].strip().split()
        meta['Ntot_cliped'] = np.float(lines_40_split[0])
        meta['clip_method'] = lines_40_split[1]
        meta['Nglobal_steps'] = np.int(lines[41].split('[')[0].strip())
        meta['N_chains'] = np.int(lines[42].split('[')[0].strip())
        meta['NEX0s_base'] = np.int(lines[43].split('[')[0].strip())
        # Synthesis Results - Best model
        meta['chi2_Nl_eff'] = np.float(lines[49].split('[')[0].strip())
        meta['adev'] = np.float(lines[50].split('[')[0].strip())
        meta['sum_of_x'] = np.float(lines[52].split('[')[0].strip())
        meta['Flux_tot'] = np.float(lines[53].split('[')[0].strip())
        meta['Mini_tot'] = np.float(lines[54].split('[')[0].strip())
        meta['Mcor_tot'] = np.float(lines[55].split('[')[0].strip())
        meta['v0_min'] = np.float(lines[57].split('[')[0].strip())
        meta['vd_min'] = np.float(lines[58].split('[')[0].strip())
        meta['AV_min'] = np.float(lines[59].split('[')[0].strip())
        meta['YAV_min'] = np.float(lines[60].split('[')[0].strip())
    except Exception:
        return meta, state
    state = True
    return meta, state


def _test_read_starlight_output_header():
    f = open('/pool/projects/starlight/STARLIGHTv04/0414.51901.393.cxt.sc4.C99.im.CCM.BN_11')
    lines = f.readlines()
    f.close()
    meta, state = read_starlight_output_header(lines)
    print(state)
    print(meta)


def _test_starlight_output():
    filepath = ('/pool/projects/starlight/STARLIGHTv04/'
                '0414.51901.393.cxt.sc4.C99.im.CCM.BN_11')
    so = StarlightOutput(filepath)
    so.pprint()
    so.write_fits('/pool/projects/starlight/STARLIGHTv04/'
                  '0414.51901.393.cxt.sc4.C99.im.CCM.BN_11.fits',
                  clobber=True)


if __name__ =='__main__':
    # _test_starlight_grid()
    # _test_starlight_base()
    # _test_starlight_config()
    # _test_starlight_mask()
    # _test_read_starlight_output_header()
    # _test_starlight_output()
    _test_starlight_mask()