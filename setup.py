import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='bopy',
    version='0.4.0',
    author='Bo Zhang',
    author_email='bozhang@nao.cas.cn',
    description='Bo Zhang (@NAOC)''s python package.', # short description
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='http://github.com/hypergravity/bopy',
    packages=setuptools.find_packages(),
    #packages=['song', 'twodspec'],
    license='New BSD',
    classifiers=["Development Status :: 5 - Production/Stable",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: MIT License",
                 "Operating System :: OS Independent",
                 "Programming Language :: Python :: 3.7",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Astronomy"],
    # package_dir={'song': 'song',
    #              'twodspec': 'twodspec'}, commented because it's wrong
    package_data={"bopy": ['bopy/data/test_spectra/*/*.fits', ],
                  "": ["LICENSE"]
                  },
    # include_package_data=False, # commented to include data!
    requires=['numpy', 'scipy', 'matplotlib', 'astropy', 'emcee', 'joblib']
)
