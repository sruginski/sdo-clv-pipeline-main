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
        'matplotlib==3.5.3',
        'sunpy==4.0.5',
        'astropy==5.2',
        'scipy==1.9.2',
        'reproject==0.9',
        'scikit-image==0.19.3',
        'pyshtools==4.10.3'
    ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='==3.9.*',
)

