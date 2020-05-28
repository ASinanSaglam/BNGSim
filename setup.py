import setuptools
try: 
    import pkg_utils
except ImportError:
    import subprocess
    import sys
    subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "pkg_utils"])
    import pkg_utils
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

name = "BNGSim"
dirname = os.path.dirname(__file__)
meta = pkg_utils.get_package_metadata(dirname, name)

setuptools.setup(
        name="BNGSim",
        version="0.3.23",
        author="Ali Sinan Saglam",
        author_email="asinansaglam@gmail.com",
        description="A simple python front-end for BioNetGen simulations",
        long_description=long_description,
        long_description_content_type="test/markdown",
        url="https://github.com/ASinanSaglam/BNGSim",
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "Operating System :: OS :: OS Independent"],
        python_requires=">=3.6",
        install_requires=meta.install_requires,
)
