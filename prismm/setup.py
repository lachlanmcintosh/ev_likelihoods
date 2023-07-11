from setuptools import setup, find_packages

setup(
    name='clonal_trees',
    version='0.1',
    description='A package to simulate cancer genomes',
    author='Your Name',
    author_email='your.email@example.com',
    packages=find_packages(),
    install_requires=[
        # list of dependencies
    ],
    entry_points={
        'console_scripts': [
            'rs=clonal_trees.run_simulation.run_simulation:main',
            'rsio=clonal_trees.run_simulation.run_simulation_IO:main',
            'do_all2=clonal_trees.clonal_trees.do_all2:main',
        ],
    },
)
