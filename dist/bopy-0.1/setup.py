from distutils.core import setup

setup(
    name = 'bopy', #*
    version = '0.1', #*
    author = 'Bo Zhang', #*
    author_email = 'bozhang@nao.cas.cn', #*
#    py_modules = ['bopy','spec','core'],
    description = 'bopy is coming!', # short description
    license = 'New BSD',
    install_requires = ['numpy>=1.7','scipy','matplotlib','nose'],
    url = 'http://github.com/hypergravity/bopy',
    classifiers = [
        "Development Status :: 6 - Mature",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Programming Language :: C",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics"],
    package_dir = {'bopy/': ''},
    packages = ['bopy','bopy/spec','bopy/core'],
    package_data = {'bopy/data':[''],
                    "": ["README.rst","README.pdf","README.dev","LICENSE","AUTHORS.rst","HISTORY.txt"]},
    include_package_data = True
    )

