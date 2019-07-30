test/bin/run-linters.sh || true
nosetests test/*.py
ngs_doit_single foo
