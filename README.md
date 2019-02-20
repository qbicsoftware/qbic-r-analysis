# R container library for reporanalysis
[![Build Status](https://travis-ci.org/qbicsoftware/r-container-lib.svg?branch=master)](https://travis-ci.org/qbicsoftware/r-container-lib)

_A collection of project-related Dockerfiles for a controlled R environment with defined R packages._

__THIS IS WORK IN PROGRESS !__

## Build new R analysis projects from the template

If you have a new R analysis project and want to add it to the R-container-lib, please follow the process described in [Rmageddon](https://github.com/qbicsoftware/r-lint-cli/blob/master/README.rst). This README assumes that you have the setup as explained in [Rmageddon](https://github.com/qbicsoftware/r-lint-cli/blob/master/README.rst).

### 1. Fork this repo

Please check the [GitHub help pages](https://help.github.com/articles/fork-a-repo/) for that.

### 2. Clone the fork

As easy as:

```
git clone https://github.com/<yourname>/r-container-lib

```

### 3. Create a new R analysis project using Rmageddon
This step is only necessary if you have not yet created a project using Rmageddon!    
Please refer to [Rmageddon](https://github.com/qbicsoftware/r-lint-cli/blob/master/README.rst) .

### 4. Modify the template

Next, you still need to make some adjustments. The current project structure looks like this:

```
projects/
    projectA/
        Dockerfile
        rpackages.txt
        scripts/
          myscript.R
    projectB/
        ...
```

When you have created your project, you have to replace the content of the `rpackages.txt` with the packages you want to have installed in your R container.

__Export the packages from the R session variable__

RCOMMAND


## Author

This repo was created by Sven Fillinger ([@sven1103](https://github.com/sven1103)), [Quantitative Biology Center](http://qbic.life), University of Tübingen.
