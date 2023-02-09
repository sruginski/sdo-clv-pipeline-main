import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sdo-pypline",
    version="0.0.1",
    author="Michael L. Palumbo",
    author_email="palumbo@psu.edu",
    description="Reduce SDO images.",
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib==3.5.1',
        'PyQt5==5.15.6',
        'pyshtools==4.10.1',
        'sunpy==4.0.4',
        'astropy==5.1',
        'scipy==1.9.2',
        'reproject==0.8',
        'scikit-image==0.19.3'
    ],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

