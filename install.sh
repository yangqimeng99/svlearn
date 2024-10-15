#!/bin/bash

PYTHON_PATH=$(realpath $(which python3))

if [ -z "$PYTHON_PATH" ]; then
    echo "Python3 is not installed or not found in PATH."
    exit 1
fi

sed -i "1s|#!.*python[0-9.]*|#!$PYTHON_PATH|" ./bin/svlearn

find ./bin/ -type f -exec chmod 755 {} \;

echo "SVLearn has been Successfully installed!"
