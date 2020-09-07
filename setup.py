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
    long_description=readme,
    author='Wolfgang Heidrich',
    author_email="wolfgang.heidrich@kaust.edu.sa",
    license="Creative Commons Atribute Non-Commercial (CC BY-NC)",
    install_requires=[
        'jupyter_client', 'IPython', 'ipykernel', 'ply'
    ],
    entry_points={
        "pygments.lexers" : ["dalexer = dacalc.pygment_lexer:DALexer"],
        "console_scripts" : ["dacalc = dacalc.calculator:main"]
    },
    classifiers=[
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3',
    ],
)


