from setuptools import setup, find_packages

setup(
    name='operon_predictor',
    version='0.0.1',
    description=(
        'predict operon'
    ),
    author='<SeraphZ>',
    author_email='<zyc199685@gmail.com>',
    license='MIT',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.14',
        'scipy',
        'matplotlib',
    ],
    url='https://github.com/Seraph-YCZhang/operon-predict',

)