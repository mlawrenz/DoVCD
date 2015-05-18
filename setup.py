from setuptools import setup
# you may need setuptools instead of distutils

setup(
    # basic stuff here
      name='DoVCD',
      author='Morgan Lawrenz',
      author_email='mlawrenz@amgen.com',
      version = '1.0',
      description = 'Process gaussian files',
    scripts = [
        'scripts/DoVCD.py'
    ]
)
