import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Metalign",
    version="0.11.0",
    author="Nathan LaPierre",
    author_email="nathanl2012@gmail.com",
    description="Metalign: efficient alignment-based metagenomic profiling via containment min hash",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nlapier2/Metalign",
    packages=['Metalign'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.5",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.5',
)
