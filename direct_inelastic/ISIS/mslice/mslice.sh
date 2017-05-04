#!/bin/sh
#
# Basic wrapper script for running in development mode. It assumes the current
# working directory is the directory containing this script.
#
PY_VERS=2.7

# Get the directory this script resides in. This will only work on bash-like shells
SCRIPT=$(readlink -f "$0")
SCRIPT_PATH=$(dirname "$SCRIPT")

INSTALL_PREFIX=/opt
PACKAGE=mantidnightly
MANTIDPYTHON=mantidpython
MANTIDPYTHON_ARGS="--classic"
MAIN_SCRIPT=${SCRIPT_PATH}/start_mslice.py

PYTHONPATH=${SCRIPT_PATH}:$PYTHONPATH ${INSTALL_PREFIX}/${PACKAGE}/bin/${MANTIDPYTHON} ${MANTIDPYTHON_ARGS} ${MAIN_SCRIPT}
