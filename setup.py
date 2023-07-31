from setuptools import setup, find_packages

setup(
    name='myopic_mces',
    version='1.0.0',
    url='https://github.com/motoharu-yano/myopic-mces',
    readme='README.md',
    license='LICENSE',
    description='A package for computation of the myopic MCES distance',
    long_description='A package for computation of the myopic MCES distance',
    platforms=['any'],
    python_requires='>=3.8',
    package_dir={"": "."},
    packages=find_packages(),
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
