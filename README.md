# CircularMapper
A method to improve mappings on circular genomes, using the BWA mapper.

[![Build Status](https://travis-ci.org/apeltzer/CircularMapper.svg?branch=master)](https://travis-ci.org/apeltzer/CircularMapper)
[![codecov](https://codecov.io/gh/apeltzer/CircularMapper/branch/master/graph/badge.svg)](https://codecov.io/gh/apeltzer/CircularMapper)

Author: Alexander Peltzer <alexander.peltzer@uni-tuebingen.de>

## Releases

The releases for this project can be found under `Releases`.

## Build the tool manually

We are using Gradle (https://gradle.org/) to build our Java applications, so you'd need to install Gradle (2+) and Java (8+) on your machine to build the executable. Now run:

```bash
git clone https://github.com/apeltzer/CircularMapper
cd CircularMapper
gradle build 
```
After that, you should find the appropriate JAR files in the subdirectories: `toolname/build/libs/`

A use-case to use this outside of the EAGER pipeline itself, is available [here](http://circularmapper.readthedocs.io/en/latest/). 
