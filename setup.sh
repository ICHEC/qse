#!/bin/bash

echo "Creating soft links for the hooks"

ln -fs -v hooks/* .git/hooks
