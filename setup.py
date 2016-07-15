from distutils.core import setup


if __name__ == '__main__':
    setup(
        name='bopy',
        version='0.2.0',
        author='Bo Zhang',
        author_email='bozhang@nao.cas.cn',
        # py_modules=['bopy','spec','core'],
        description='Bo Zhang (@NAOC)''s python package.',  # short description
        license='New BSD',
        # install_requires=['numpy>=1.7','scipy','matplotlib','nose'],
        url='http://github.com/hypergravity/bopy',
        classifiers=[
            "Development Status :: 6 - Mature",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Operating System :: OS Independent",
            "Programming Language :: C",
            "Programming Language :: Python :: 2.7",
            "Topic :: Scientific/Engineering :: Astronomy",
            "Topic :: Scientific/Engineering :: Physics"],
        package_dir={'bopy/': ''},
        packages=['bopy',
                  'bopy/spec',
                  'bopy/core',
                  'bopy/mcmctools',
                  'bopy/obstools',
                  'bopy/obstools/Xinglong216HRS',
                  'bopy/helpers',
                  'bopy/helpers/ezpadova'
                  ],
        package_data={'bopy/data': [''],
                      "":          ["LICENSE"]},
        include_package_data=True,
        requires=['numpy', 'scipy', 'matplotlib', 'astropy', 'lmfit']
    )
