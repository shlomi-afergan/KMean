from setuptools import setup, find_packages, Extension



# setup() parameters - https://packaging.python.org/guides/distributing-packages-using-setuptools/
setup(
    name='mykmeansspf',
    version='0.1.0',
    author="Ido and Shlomi",
    author_email="shlomido@yesmail.com",
    description="kmeanssp C-API",
    install_requires=['invoke'],
    packages=find_packages(),  # find_packages(where='.', exclude=())
    #    Return a list of all Python packages found within directory 'where'
    license='GPL-2',
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3 :: Only',
        # We need to tell the world this is a CPython extension
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            # the qualified name of the extension module to build
            'mykmeansspf',
            # the files to compile into our module relative to ``setup.py``
            sources=['spkmeansmodule.c','spkmeans.c']
        )
    ]
)
