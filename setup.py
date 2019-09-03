import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="xyz2graph",
    version="0.1",
    author="Mykola Zotko",
    description="Package for reading of .xyz files and constructing of molecular graphs from atomic coordinates.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zotko/xyz2graph.py",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['plotly', 'networkx'],
    python_requires='>=3.5',
)