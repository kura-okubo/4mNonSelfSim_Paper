#!bin/sh
# https://stackoverflow.com/a/5200267
find . -size +10M | cat >> .gitignore
