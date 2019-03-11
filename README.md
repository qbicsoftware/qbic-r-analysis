# R container library for reporanalysis
[![Build Status](https://travis-ci.org/qbicsoftware/r-container-lib.svg?branch=master)](https://travis-ci.org/qbicsoftware/r-container-lib)

_A collection of project-related Dockerfiles for a controlled R environment with defined R packages._

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
        scripts/
          myscript.R
        data/
        Dockerfile
        environment.yml     
    projectB/
        ...
```

Add a subfolder project[X] for your new project and copy its contents into it. Ensure that it contains a Dockerfile, an environment.yml file *after* the building process of Rmageddon, a data folder and finally a scripts folder with your R scripts.

### 5. Sanity check: lint

If you want to finally verify that your newly added project doesn't break any requirements you can run Rmageddon lint again on it. Please refer to [Rmageddon](https://github.com/qbicsoftware/rmageddon-cli/blob/master/doc/Rmageddon.md) .

### 6. Publish your R project

The only thing left is to submit a pull request or directly commit to this repository. 
If you don't remember git so well, this condensed cheatsheet may help [Git Cheatsheet](https://www.keycdn.com/blog/git-cheat-sheet).

## Author

This repo was created by Sven Fillinger ([@sven1103](https://github.com/sven1103)), [Quantitative Biology Center](http://qbic.life), University of Tübingen.
