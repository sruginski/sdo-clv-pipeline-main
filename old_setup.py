import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sdo_clv_pipeline",
    version="0.1.0",
    author="Michael L. Palumbo",
    author_email="palumbo@psu.edu",
    description="SDO data reduction for center-to-limb variability",
    install_requires=[
        'numpy',
        'pandas',
        'PyQt5',
        'matplotlib',
        'sunpy',
        'astropy',
        'scipy',
        'reproject',
        'scikit-image',
        'pyshtools'
    ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
)

