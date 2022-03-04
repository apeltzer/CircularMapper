# CircularMapper
A method to improve mappings on circular genomes, using the BWA mapper.

[![codecov](https://codecov.io/gh/apeltzer/CircularMapper/branch/master/graph/badge.svg)](https://codecov.io/gh/apeltzer/CircularMapper)

:warning: CM is only tested with `bwa aln` and _not_ with the `bwa mem` method. You might experience difficulties due to differences how `bwa mem` handles CIGAR strings in SAM/BAM files and thus are advised to only use `bwa aln` mapped output. 
:information_source: As this was part of my PhD (> 5yrs ago), I'm not actively maintaining this anymore. If you find bugs, you may open issues but I will not guarantee any bugfixes or even looking at them anymore - code is open source, so if you can and want to fix things, I am however happy to have a look at a PR to fix things. 

Author: Alexander Peltzer

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
