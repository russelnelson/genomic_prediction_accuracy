#!/bin/bash

#finds files that are larger than 50mb and adds them to the .gitignore

find . -size +50M -exec ls {} \+ | grep -v ".git/" - | cat - .gitignore | sort | uniq > out
rm .gitignore
mv out .gitignore
