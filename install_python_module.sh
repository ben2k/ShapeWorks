#!/bin/bash -x

# Installs shapeworks Python module
# Ensures conda Python module is able to locate correct shapeworks_py binary module.

# install pip module
if ! pip install Python/shapeworks; then exit -1; fi

CONDA_INSTALL_DIR=`pip show shapeworks | grep Location | awk '{print $2}'`/shapeworks

python ./update_installed_python_module_binary_path.py $CONDA_INSTALL_DIR/setup_path.py `pwd`/bin
