from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Neoantigens from Structural Variations'
LONG_DESCRIPTION = 'neosv is a workflow for annotating SVs and predicting neoantigens from SVs'
# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="neosv", 
        version=VERSION,
        author="Yang Shi",
        author_email="<shiyang_bio@163.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=['pyensembl>=1.9.3', 'biopython>=1.79'],
        entry_points={
            "console_scripts": [
                'neosv=neosv.main:main'
            ]
        },
        keywords=['python', 'bioinformatics', 'genomics', 'neoantigen', 'structural variation'],
        classifiers= [
            "Development Status :: 5 - Production/Stable",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS",
            "Operating System :: OS/2"
        ]
)
