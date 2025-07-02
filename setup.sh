#!/bin/bash

echo "Copying the hooks to .git/hooks/"

for x in `ls -1 hooks/*`
do
    cp -f -v -p ./$x ./.git/$x
done
