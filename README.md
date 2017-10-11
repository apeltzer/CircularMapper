# CircularMapper
A method to improve mappings on circular genomes, using the BWA mapper.

[![Build Status](https://lambda.informatik.uni-tuebingen.de/jenkins/buildStatus/icon?job=CircularMapper)](https://lambda.informatik.uni-tuebingen.de/jenkins/view/EAGER/job/CircularMapper/)
[![codecov](https://codecov.io/gh/apeltzer/CircularMapper/branch/master/graph/badge.svg)](https://codecov.io/gh/apeltzer/CircularMapper)
[ ![Download Realigner](https://api.bintray.com/packages/apeltzer/EAGER/CircularMapper.Realigner/images/download.svg) ](https://bintray.com/apeltzer/EAGER/CircularMapper.Realigner/_latestVersion)
[ ![Download Generator](https://api.bintray.com/packages/apeltzer/EAGER/CircularMapper.Generator/images/download.svg) ](https://bintray.com/apeltzer/EAGER/CircularMapper.Generator/_latestVersion)

Author: Alexander Peltzer <alexander.peltzer@uni-tuebingen.de>

## Releases

The releases for this project can be found on [Bintray](https://bintray.com/apeltzer/EAGER/) or direct download from there [respectively](https://dl.bintray.com/apeltzer/EAGER/com/uni-tuebingen/de/it/eager/).

## Build the tool manually

We are using Gradle (https://gradle.org/) to build our Java applications, so you'd need to install Gradle (2+) and Java (8+) on your machine to build the executable. Open the ``build.gradle`` file and delete both `bintray { ... }` sections (once for each tool, circulargenerator and realignsamfile) we use for deploying the CI build to Bintray normally. Save the file and next you can run:

```bash
git clone https://github.com/apeltzer/CircularMapper
cd CircularMapper
gradle build 
```
After that, you should find an executable JAR file in the directories `generator/build/libs` and `realigner/build/libs`. 

A use-case to use this outside of the EAGER pipeline itself, is available [here](http://circularmapper.readthedocs.io/en/latest/). 
