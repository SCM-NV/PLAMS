# Contributing to PLAMS

Thank you for your interest in contributing to the PLAMS project! 

We want to keep it as easy as possible to contribute changes that
get things working in your environment. There are a few guidelines that we
need contributors to follow so that we can have a chance of keeping on
top of things.

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free)
* [Fork](https://github.com/SCM-NV/PLAMS/fork) the repository on GitHub

## Making Changes

* Create a feature branch from where you want to base your work
  * This is usually the trunk branch - only target release branches if you are certain your fix must be on that
    branch
  * To quickly create a feature branch based on trunk run `git checkout -b
    MyName/my_contribution trunk`
* Make commits of logical units
* Make sure your commit messages are informative
* Make sure you have added the necessary tests for your changes

## Submitting Changes

* Push your changes to a feature branch in your fork of the repository
* Submit a pull request to the repository in the SCM-NV organization
* The core team reviews Pull Requests on a regular basis and will provide feedback / approval

## Developer Tools

In order to maintain a high code quality, we use a variety of code developer tools in our repo.
These include:

| Tool                                        | Purpose                                  |
|:--------------------------------------------|:-----------------------------------------|
| [black](https://black.readthedocs.io/)      | code formatting                          |
| [mypy](https://mypy.readthedocs.io/)        | static type checking                     |
| [pytest](https://docs.pytest.org/)          | testing                                  |
| [coverage](https://coverage.readthedocs.io) | measuring test coverage                  |
| [sphinx](https://www.sphinx-doc.org/)       | generating documentation                 |

These can be installed with the command `pip install '.[dev]'`.

A series of checks are conducted in our CI pipeline, which must pass 
before your pull request is considered ready for review.