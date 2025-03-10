from setuptools import setup, find_packages

# Check if libbdsg is already installed
try:
    import bdsg
except ImportError:
   raise ImportError("ERROR bdsg not found")

setup(
    name="STOAT",
    version="0.1.0",
    author="Matis Alias-Bagarre",
    author_email="matis84700@gmail.com",
    description="STOAT is a specialized tool developed for conducting Genome-Wide Association Studies (GWAS) with a unique focus on snarl structures within pangenome graphs.",
    long_description=open("README.md").read(),
    url="https://github.com/Plogeur/STOAT",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=[
        "cyvcf2==0.31.1",
        "numpy==2.0.1",
        "pandas==2.2.2",
        "statsmodels==0.14.4",
        "qmplot==0.3.3",
        "scipy==1.14.1",
        "setuptools==75.6.0",
        "scikit-learn==1.6.0",
        "matplotlib==3.9.4",
        "seaborn==0.13.2",
        "plotly==5.24.1",
        "pytest==8.3.4"
    ],
    entry_points={
        "console_scripts": [
            "stoat=stoat:main",  # CLI entry point
        ]
    },
)
