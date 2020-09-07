from setuptools import setup, find_packages
import sys
from os import chdir

with open('README.md') as f:
    readme = f.read()


print('-'*50,file=sys.stderr)
print("Packages:",find_packages(),file=sys.stderr)
print('-'*50,file=sys.stderr)

setup(
    name='dacalc',
    version='1.0',
    packages=find_packages(),
    description='Dimensional Analysis Calculator',
    author='Wolfgang Heidrich',
    author_email="wolfgang.heidrich@kaust.edu.sa",
    license="Creative Commons Atribute Non-Commercial (CC BY-NC)",
    url="https://github.com/wgheidrich/DACalc",
    install_requires=[
        'jupyter_client', 'IPython', 'ipykernel', 'ply'
    ],
    entry_points={
        "pygments.lexers" : ["dalexer = dacalc.pygment_lexer:DALexer"],
        "console_scripts" : ["dacalc = dacalc.calculator:main"]
    },
    classifiers=[
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Manufacturing',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Programming Language :: Python :: 3',
    ],
    long_description=readme,
)


