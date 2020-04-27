import os
import setuptools

SCRIPTS = []
SCRIPTS.extend([os.path.join("scripts", script)
				for script in os.listdir(os.path.join(os.path.dirname(__file__), "scripts"))
				if script.endswith(".py")])

HERE = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(HERE, 'README.md'), 'r') as fid:
	LONG_DESCRIPTION = fid.read()

setuptools.setup(
    name="Metalign",
    version="0.12.1",
    author="Nathan LaPierre",
    author_email="nathanl2012@gmail.com",
    description="Metalign: efficient alignment-based metagenomic profiling via containment min hash",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/nlapier2/Metalign",
    #packages=['Metalign'],
    package_data={'Metalign': ['data/cmash*', 'data/db_info*', 'data/organism_files/*']},
    scripts=SCRIPTS,
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
