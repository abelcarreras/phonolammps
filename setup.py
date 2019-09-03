try:
    from setuptools import setup, Extension
    use_setuptools = True
    print("setuptools is used.")
except ImportError:
    from distutils.core import setup, Extension
    use_setuptools = False
    print("distutils is used.")


def get_version_number():
    for l in open('phonolammps/__init__.py', 'r').readlines():
        if not(l.find('__version__')):
            exec(l, globals())
            return __version__

setup(name='phonoLAMMPS',
      version=get_version_number(),
      description='phonoLAMMPS module',
      author='Abel Carreras',
      url='https://github.com/abelcarreras/phonolammps',
      author_email='abelcarreras83@gmail.com',
      packages=['phonolammps'],
      scripts=['scripts/phonolammps'],
      install_requires=['phonopy', 'numpy', 'matplotlib', 'seekpath', 'dynaphopy'],
      license='MIT License'
      )
