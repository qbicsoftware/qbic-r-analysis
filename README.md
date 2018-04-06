# QBiC R analysis
_A collection of project-related Dockerfiles for a controlled R environment with defined R packages._

__!! THIS IS WORK IN PROGRESS !! __

## Build new R analysis projects from the template

We have a R analysis template and you can use the `cookiecutter` tool for that. Please check out [https://github.com/qbicsoftware/qbic-r-analysis-template](https://github.com/qbicsoftware/qbic-r-analysis-template) and follow the instructions shown in the README.

### Fork this repo

Please check the [GitHub help pages](https://help.github.com/articles/fork-a-repo/) for that.

### Clone the fork

As easy as:

```
git clone https://github.com/<yourname>/qbic-r-analysis
```

### Create a new R analysis project using cookiecutter

Use cookiecutter and a project from a template with:

```
cookiekutter https://github.com/qbicsoftware/qbic-r-analysis-template
```

You will get asked to answer some questions, and cookiecutter will automatically put the information into the template!

```
> cookiekutter https://github.com/qbicsoftware/qbic-r-analysis-template
r_version [3.2.4]: 
author_name [Sven Fillinger]: 
author_email [sven.fillinger@qbic.uni-tuebingen.de]: 
container_version [0.1dev]: 
project_code [QABCD]: 
```
The values in `[]` are the __default values__, that are taken if you just hit `enter`.

### Modify the template

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
