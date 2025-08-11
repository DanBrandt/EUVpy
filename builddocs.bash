#!/bin/bash
# Installs and builds Sphinx with necessary extensions. If run more than once, will search for updates in codebase to add to html
# 05/2022 Aidan Sorensen

#venvName="venv"
#if [[ -d $venvName ]]; then
#    echo "Virtual environment appears to already exist."
#else
#    echo "Creating new virtual environment"
#    python3.7 -m venv $venvName
#fi
#
##it does need to be activated, right?
#source venv/bin/activate

# pip3 install -r ./sphinx_requirements.txt
# variables
author="Daniel Brandt"
release=1
version=1.0.0
currentdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )" # Find current directory
echo "Current dir is: $currentdir"

mkdir -p docs/source
cp -r $currentdir/staticdocs docs/source

# remove pre-existing docs before rebuilding
rm -rf docs/source/*.rst
rm -rf docs/source/conf.py
rm -rf docs/build/*

OLDPYTHONPATH=$PYTHONPATH

PYTHONPATH=$currentdir/src/:$PYTHONPATH #add project dir to python path
export PYTHONPATH
echo "The PYTHONPATH is:"
echo $PYTHONPATH
proj_name=$(basename $currentdir) # Project name generated from the parent folder name

# subshell
(
    if [ ! -d docs/ ]; then
        mkdir -p docs/
    fi
    cd docs/

#    FILE=source/conf.py
#    if test -f "$FILE"; then
#        echo "$FILE already exists, skipping sphinx-quickstart."
#    else
#        # sphinx-quickstart is a tool that asks questions about the project, then generates documentation dir and makefile for sphinx-build
#        sphinx-quickstart \
#          -q `# -q supresses extra config settings` \
#          --sep `# --sep separates build and source directory inside docs` \
#          -a="$author" `# -a Author name var 'author' set at top of script` \
#          -p=$proj_name `# -p Project name from parent dir` \
#          -v=$version `# -v Version of the project, set to 1` \
#          -r=$release `# -r Release number, set to 1` \
#          --language=en `# --language Sets language, en="English"` \
#          --ext-autodoc `# --ext-autodoc Extension: For use with apidoc command, employs autodoc to search pythonic package structure for docstrings` \
#          --ext-mathjax `# --ext-mathjax Extension: for use specifically with LaTeX formatting, mostly unnecessary given .. math:: function` \
#          --ext-viewcode `# --ext-viewcode Extension: Enables source code to be viewed from generated html` \
#          --makefile `# --makefile enables the generation of a makefile` \
#          --extensions sphinx.ext.napoleon


    #     # Hacky solution to modify the contents of conf.py. Should look for a nicer Bashy solution
    #     echo -e ''"import os\nimport sys\nsys.path.insert(0, os.path.abspath('../src/'))\n"'' > path_extension.txt
    #     cat path_extension.txt $FILE  > source/updated_conf.py
    #     mv source/updated_conf.py source/conf.py





    # #     # Add Table of Contents to index.rst
    # #     # ${proj_name}
    # #     # touch tmp
    # #     # awk -v n=9 -v s=".. toctree::\n  :maxdepth: 4\n  :caption: Contents:\n\n  modules\n\n" \
    # #     #                 "NR == n {print s} {print}" \
    # #     #                 source/index.rst > source/tmp; mv source/tmp source/index.rst
#    fi

    # # sphinx-apidoc is a command that uses autodoc to import/find docstrings in all packages/subpackages from parent dir
    sphinx-apidoc \
    --templatedir=templates/apidoc \
     -F `# -F forces the generation of conf.py` \
     -a `# -a append module path to sys path` \
     -e `# -e generates a separate rst file for each module, so they show up on separate pages` \
     -M `# -M shows modules first, then submodules. (less redundancy in the TOC)` \
     -A "$author" `# -A Author name var 'author' set at top of script` \
     -H $proj_name `# -H Project name from parent dir` \
     -V $version `# -V Version of the project, set to 1` \
     -R $release `# -R Release number, set to 1` \
     -o source/ `# -o source/ tells apidoc to output to the given path (source/)` \
      ../src/ \
      --extensions sphinx.ext.napoleon,sphinx.ext.mathjax,myst_parser \



    # # Now that the conf.py file has been generated, we'll modify it to do what we want it to do:


    # sphinx-build is a command used to build the .rst files into html, put them in the build dir
    # -D modifies conf.py at runtime to include these options. Does not save these changes in conf.py for later use
    sphinx-build -b html `# -b build in this format:html` \
     -D html_theme="sphinx_rtd_theme" `# Sets theme to readthedocs, looks significantly better than alabaster` \
     -D html_theme_options.sticky_navigation=True `# Pins navigation so it stays on page when scrolling` \
     -D html_theme_options.collapse_navigation=False `# Allows expandable branch button in TOC when navigating ` \
     `#-D numfig=False -D math_numfig=False -D math_number_all=False # numfig, math_numfig, math_number_all work together to enable equation labeling` \
     source/ build/html/

)

#to open static site after build: firefox docs/build/html/index.html (any .html in this dir will work)
# -D sys.path.insert(0, os.path.abspath('.')) `#add system path ?` \

PYTHONPATH=$OLDPYTHONPATH
export PYTHONPATH
