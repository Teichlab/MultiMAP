from setuptools import setup, find_packages

setup(
	name='MultiMAP',
	version='0.0.1',
	description='MultiMAP',
	url='https://github.com/Teichlab/MultiMAP',
	packages=find_packages(exclude=['docs', 'examples']),
	install_requires=['numpy','scipy','numba','scikit-learn'],
	author='Mika Sarkin Jain',
	author_email='mikasarkinjain@gmail.com',
	license='MIT'
)
