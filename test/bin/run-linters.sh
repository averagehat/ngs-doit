#!/bin/bash

mypy ngs_doit/*.py
mypy test/*.py
mypy ./*.py

pylint ngs_doit test

test/bin/shellcheck *.sh
test/bin/shellcheck test/bin/*.sh
