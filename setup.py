from setuptools import setup, find_packages

setup(
    name="dbfold2",
    version="0.1.0",
    author="Kibum Park",
    author_email="your-email@example.com",
    description="A project description here.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/dbfold2",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    # install_requires=[],  # Omit if you don't want to specify dependencies
)