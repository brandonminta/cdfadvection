from setuptools import setup, find_packages

setup(name='cfdadvection',
      description='CDS resolutions with different slope limiters',
      author='Brandon Minta',
      author_email='brandon.minta@yachaytech.edu.ec,
      license='MIT',
      packages=find_packages(),
      install_requires=['numpy', 'matplotlib'])