from setuptools import setup

setup(
    name="relatively",
    version="1.0.4",
    description="Abundance figure with absolute and relative abundances across a hierarchy",
    url="https://github.com/brwnj/relatively",
    author="Joe Brown",
    author_email="brown@pm.me",
    license="MIT",
    packages=["relatively"],
    install_requires=["numpy", "pandas>=0.23.1", "plotly>=3.0.0"],
    zip_safe=False,
)
