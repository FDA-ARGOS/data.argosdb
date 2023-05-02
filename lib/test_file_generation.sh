#!/bin/bash

# Simple script to take the HEAD for test files.
# Need to provide 2 directories as the arguments.
# first dir is the location of files to test
# second dir is the location fo store the test files.

for i in $(ls $1)
do
    head $1$i > $2$i
done