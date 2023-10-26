from setuptools import setup, find_packages

with open('requirements.txt', 'r') as f:
    requirements = f.read().splitlines()

setup(name='metamotif',
      version='0.5.0',
      description='metamotif',
      url='http://github.com/mhorlacher/metamotif',
      author='Marc Horlacher',
      author_email='marc.horlacher@helmholtz-muenchen.de',
      license='MIT',
      install_requires=requirements,
      packages=find_packages(),
      include_package_data=True,
      entry_points = {
            'console_scripts': [
                  'metamotif=metamotif.__main__:main',
            ],
      },
      zip_safe=False)
