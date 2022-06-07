from setuptools import setup, find_packages

requirements = [
    "pandas",
    "numpy",
    "logomaker",
    "sklearn",
    "umap-learn",
]

setup(name='metamotif',
      version='0.1',
      description='metamotif',
      url='http://github.com/mhorlacher/metamotif',
      author='Marc Horlacher',
      author_email='marc.horlacher@helmholtz-muenchen.de',
      license='MIT',
      install_requires=requirements,
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False)
