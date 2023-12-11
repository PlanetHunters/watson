import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = "0.8.1"
setuptools.setup(
    name="dearwatson", # Replace with your own username
    version=version,
    author="M. Dévora-Pajares",
    author_email="mdevorapajares@protonmail.com",
    description="Visual Vetting and Analysis of Transits from Space ObservatioNs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/PlanetHunters/watson",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=[
                        "bokeh==2.4.2", # TPFPlotter dependency
                        'configparser==5.0.1',
                        "extension-helpers==0.1",
                        "exoml==0.0.15",
                        "imageio==2.9.0",
                        "lcbuilder==0.14.2",
                        "matplotlib==3.8.2",
                        'pyparsing==2.4.7', # Matplotlib dependency
                        "pyyaml==6.0.1",
                        "pillow==9.5.0",
                        "reportlab==4.0.7",
                        'setuptools>=41.0.0'
    ]
)
