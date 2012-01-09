#!/bin/bash

PYTHONPATH=$PYTHONPATH:"$(pwd)/src"

python2 -m unittest discover -s tests

