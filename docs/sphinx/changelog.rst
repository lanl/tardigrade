.. _changelog:


#########
Changelog
#########


******************
0.1.9 (unreleased)
******************

Internal Changes
================
- Added SSL workaround to Pages job (:issue:`23`, :merge:`29`). By `Sergio Cordova`_.

******************
0.1.8 (2022-04-28)
******************

Internal Changes
================
- Remove the package deployment Gitlab-CI job because the AEA Compute environment no longer allows projects to directly
  update the environment. Instead, projects must request that their package is added to the AEA Compute environment
  build (:issue:`18`, :merge:`23`). By `Kyle Brindley`_.
- Move the production release automatic microbumping to a dedicated Gitlab-CI job (:issue:`18`, :merge:`23`). By `Kyle
  Brindley`_.
- Added fix to avaid warnings treated as errors introduced in Sphinx 5 (:issue:`25`, :merge:`30`). By `Sergio Cordova`_.


******************
0.1.7 (2022-03-24)
******************

Internal Changes
================
- Test and deploy against the "aea-release" and "aea-beta" environments for pending AEA Compute Environment changes:
  https://ddw-confluence.lanl.gov/display/PYT/2022/02/08/AEA+Compute+environment+updates+coming+March+31%2C+2022
  (:merge:`21`). By `Kyle Brindley`_.


******************
0.1.6 (2022-03-21)
******************

Bug fixes
=========
- Update the documentation ``cmake`` command to match the new documentation directory structure (:merge:`10`). By `Kyle
  Brindley`_.
- Re-enabled the Abaqus integration tests (:merge:`14`). By `Nathan Miller`_.

Documentation
=============
- Update URLs for cpp stub repository (:issue:`22`, :merge:`28`). By `Prabhu Khalsa`_.
- Deploy both ``master`` and ``dev`` branch documentation (:issue:`4`, :merge:`8`). By `Kyle Brindley`_.
- Fix broken documentation URLs in README (:merge:`11`). By `Kyle Brindley`_.
- Fix broken Gitlab documentation URLs in Gitlab setup (:merge:`12`). By `Kyle Brindley`_.
- Fix broken ``rename`` command in Gitlab setup (:merge:`13`). By `Kyle Brindley`_.

Internal Changes
================
- Removed unused myst-parser extension from the Sphinx configuration (:issue:`9`, :merge:`15`). By `Kyle Brindley`_.
- Update the build configuration to handle conda environments than manage cpp compilers and libraries (:issue:`11`
  :merge:`16`). By `Kyle Brindley`_.
- Add back compiler flags related to code warnings for the project wide compile options (:issue:`12`, :merge:`18`). By
  `Kyle Brindley`_.


******************
0.1.5 (2021-07-19)
******************

Documentation
=============
- Update project setup instructions from Atlassian to Gitlab workflows (:issue:`2`, :merge:`4`). By `Kyle Brindley`_.

Internal Changes
================
- Convert README from markdown to restructured text (:issue:`2`, :merge:`4`). By `Kyle Brindley`_.
- Separate Abaqus integration test setup from Abaqus integration ctest declaration. Enables documentation build
  dependencies on Abaqus integration test input files without requiring Abaqus test execution on systems with no Abaqus
  installation (:issue:`2`, :merge:`4`). By `Kyle Brindley`_.


******************
0.1.4 (2021-07-13)
******************

Internal Changes
================
- Upstream project settings update to set default merge-request branch. By `Kyle Brindley`_.

******************
0.1.3 (2021-07-13)
******************

- Migrate from ddw-bibucket.lanl.gov to re-git.lanl.gov and convert to Gitlab CI/CD (:issue:`1`, :merge:`1`). By `Kyle
  Brindley`_.

******************
0.1.2 (2021-07-01)
******************

Internal Changes
================
- Use Git SCM tags for semantic versioning (:jira:`702`, :pull:`50`). By `Kyle Brindley`_.
- Master branch production release logic for CD, including automated micro-version bumps (:jira:`702`, :pull:`50`). By `Kyle
  Brindley`_.


******************
0.1.1 (2021-06-15)
******************

Bug Fixes
=========
- Corrected bug in `cpp_stub.cpp` in the map of `ddsdde` to `DDSDDE` due to using `spatialDimensions` instead
  of `NTENS` (:jira:`685`, :pull:`47`). By `Nathan Miller`_.

Documentation
=============
- Add camelCase project name replacement instructions to project setup. By `Kyle Brindley`_.


******************
0.1.0 (2021-05-28)
******************

New Features
============
- Add CMake install configuration and CI/CD scripts for build, test, and installation to a Conda environment
  (:jira:`654`, :pull:`41`). By `Kyle Brindley`_.

Documentation
=============
- Update the Python package dependencies and add an example approach to future updates to the documentation
  (:jira:`636`, :pull:`37`). By `Kyle Brindley`_.
- Add file renaming commands to the project setup instructions (:jira:`634`, :pull:`38`). By `Kyle Brindley`_.
- Update the user manual to reflect required environment variable ``LD_LIBRARY_PATH`` (:jira:`662`, :pull:`43`). By
  `Kyle Brindley`_.

Internal Changes
================
- Update markdown syntax in README for wider compatibility (:jira:`604`, :pull:`36`). By `Kyle Brindley`_.
- Maintenance on ReST style guide updates (:jira:`604`, :pull:`36`). By `Kyle Brindley`_.
- Address BOOST output test stream deprecations and update minimum version
  (:jira:`654`, :pull:`41`). By `Kyle Brindley`_.
- Change project UMAT library name to avoid conflicts with external projects (:jira:`661`, :pull:`42`). By `Kyle
  Brindley`_.
- Remove the ``CXX`` compiler variable settings for build scripts (:jira:`671`, :pull:`44`). By `Kyle Brindley`_.

Enhancements
============
- Add multi-host and multi-environment CI/CD (:jira:`630`, :pull:`39`). By `Kyle Brindley`_.


******************
0.0.4 (2021-04-30)
******************

Documentation
=============
- Clarify behavior for custom target for the integration tests (:jira:`557`, :pull:`29`). By `Kyle Brindley`_.
- Add template documentation for the Abaqus material input definition (:jira:`575`, :pull:`31`). By `Kyle Brindley`_.
- Major overhaul of documentation organization to single source the Jenkins setup information from markdown files.  Adds
  the ``myst-parser`` Python package dependency and a pull request reviewer guide (:jira:`601`, :pull:`33`). By `Kyle
  Brindley`_.

Internal Changes
================
- Update Jenkins CI configuration to build and test for PRs to both ``master`` and ``dev`` branches (:jira:`544`,
  :pull:`26`). By `Kyle Brindley`_.
- Minor cleanup to root directory files. Move configuration and environment files to a subdirectory (:jira:`544`,
  :pull:`26`). By `Kyle Brindley`_.
- Add integration test CMake target for conditional rebuilds and file copy (:jira:`551`, :pull:`27`). By `Kyle
  Brindley`_.
- Add one ctest per abaqus input file (:jira:`551`, :pull:`27`). By `Kyle Brindley`_.
- Accept paths for input file in integration test shell script and check for errors in the abaqus stdout/stderr log
  (:jira:`551`, :pull:`27`). By `Kyle Brindley`_.
- Enable parallel CMake builds for continuous integration (CI) tests (:jira:`518`, :pull:`28`). By `Kyle Brindley`_.
- Add c++ source files ``*.cpp`` as dependencies for the Doxygen CMake target (:jira:`569`, :pull:`30`). By `Kyle
  Brindley`_.
- Add checks for ``STATEV`` and ``PROPS`` vector lengths to the abaqus interface. Throw exceptions with file and
  function name to interrupt Abaqus execution on input errors (:jira:`575`, :pull:`31`). By `Kyle Brindley`_.
- Add Abaqus interface unit tests for checking the ``STATEV`` and ``PROPS`` vector lengths (:jira:`575`, :pull:`31`). By
  `Kyle Brindley`_.
- Add unit tests for error codes in ``cpp_stub::sayHello`` (:jira:`334`, :pull:`32`). By `Kyle Brindley`_.

Enhancements
============
- Add error reporting to the Abaqus interface from the ``error_tools`` package (:jira:`334`, :pull:`32`). By `Kyle Brindley`_.


******************
0.0.3 (2021-04-13)
******************

Internal Changes
================
- Use ``abaqus_tools`` from a dedicated project (:jira:`535`, :pull:`23`). By `Kyle Brindley`_.
- Add ``bibtex_bibfiles`` variable to Sphinx configuration for newer version of ``sphinxcontrib.bibtex`` extension in
  Anaconda 2020 (:jira:`526`, :pull:`21`). By `Kyle Brindley`_.
- Add explicit list of documentation source files for better conditional CMake documentation re-builds (:jira:`526`,
  :pull:`21`). By `Kyle Brindley`_.


******************
0.0.2 (2021-02-11)
******************

Breaking changes
================
- Remove testing and support for intel ``icpc`` compiler (:jira:`516`, :pull:`9`). By `Kyle Brindley`_.

New Features
============
- Add do-nothing template c++ Abaqus UMAT interface and sample Abaqus input file (:jira:`502`, :pull:`6`). By `Kyle Brindley`_.
- Use example c++ library in Abaqus UMAT template (:jira:`505`, :pull:`8`). By `Kyle Brindley`_.
- Add c++ to fortran variable conversion and Abaqus variable return template (:jira:`521`, :pull:`15`, :pull:`16`). By
  `Kyle Brindley`_.
- Add common abaqus tensor handling tools and a c++ converted umat interface (:jira:`522`, :pull:`17`). By `Kyle
  Brindley`_.

Bug fixes
=========

Documentation
=============
- Add changelog to documentation (:jira:`450`, :pull:`11`). By `Kyle Brindley`_.
- Add direct CMake build instructions and minimal user manual (:jira:`519`, :pull:`12`). By `Kyle Brindley`_.
- Add release guidance and release branch instructions (:jira:`520`, :pull:`13`). By `Kyle Brindley`_.

Internal Changes
================
- Use BOOST and ctest for unit testing (:jira:`357`, :pull:`4`). By `Kyle Brindley`_.
- Update Jenkins CI configuration and store with version controlled repository (:jira:`442`, :pull:`5`). By `Kyle Brindley`_.
- Demonstrate c++ ``vector_tools`` library for unit testing (:jira:`506`, :pull:`7`). By `Kyle Brindley`_.
- Add integration tests for Abaqus UMAT interface (:jira:`504`, :pull:`10`). By `Kyle Brindley`_.
- Move project Abaqus interface into project files. Treat UMAT Fortran/c++ subroutine as a UMAT selection and pass
  through subroutine (:jira:`523`, :pull:`18`). By `Kyle Brindley`_.
- Bump micro version number for release (:jira:`524`). By `Kyle Brindley`_.

Enhancements
============


******************
0.0.1 (2020-10-26)
******************

Breaking changes
================

New Features
============
- Create c++ stub repository targeting constitutive modeling (:jira:`332`, :pull:`1`). By `Kyle Brindley`_.

Bug fixes
=========

Documentation
=============

Internal Changes
================
- Add continuous integration scripts (:jira:`333`, :pull:`2`). By `Kyle Brindley`_.

Enhancements
============
