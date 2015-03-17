try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'name': 'overlapy',
    'description': 'A simple library to perform structural alignment on two PDB files.',
    'author': 'Alvaro Abella, Josep Arus-Pous and Samuel Miravet-Verde',
    'url': 'http://overlapy.undeadpixels.net',
    'download_url': 'https://github.com/smv818vms/overlapy',
    'author_email': 'josep.arus01@estudiant.upf.edu',
    'version': '0.0.1',
    'license': "MIT",
    'install_requires': ['Biopython', 'numpy'],
    'packages': ['overlapy'],
    'scripts': ['bin/overlapy']
}

setup(**config)
