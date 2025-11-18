""" dynamic (thermodynamic) flux balance analysis ""

.. moduleauthor:: remind team

"""

from setuptools import setup, find_packages

version_tag = '0.0.1'


setup(name='remind',
      version=version_tag,
      author='remind team',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/remind/',
      download_url='https://github.com/EPFL-LCSB/remind/archive/'+version_tag+'.tar.gz',
      #todo need to include pytables somehow or
      install_requires=['numpy',
                        'pandas',
                        'deap',
                        'optlang==1.6.0',
                        'matplotlib',
                        'numexpr >= 2.6.2',
                        'Cython >= 0.21',
                        'wheel',
                        'tables==3.6.1',
                        'tqdm',
                        # 'sklearn',
                        'seaborn',
                        'jedi == 0.17.2',  # for ipython in python 3.6
                        'cobra == 0.17.0',
                        'ipython',
                        'scikit-learn',
                        # 'bokeh==2.3.3',
                        # 'holoviews',
                        'sympy==1.3',
                        'wheel',
                        'pytfa'
                        ],
      packages = find_packages(),
      python_requires='>=3, <4',
     # python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='remind',
      keywords=['microbiome','community','models','optimal'],
      license='Apache2',

      # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry'
            'Environment :: Console',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: Apache Software License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
           # 'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
      ],
     )
