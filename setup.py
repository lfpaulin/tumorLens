from setuptools import setup, find_packages

setup(
    name='TumorLens',
    version='1.0.250801',
    packages=find_packages(),
    url='https://github.com/lfpaulin/tumorLens',
    license='MIT',
    author='Luis Paulin',
    author_email='lfpaulin@gmail.com',
    description='Somatic genetic and epigenetic analysis for cancer',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
