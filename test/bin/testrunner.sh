test/bin/run-linters.sh || true
py.test --cov=vartable test/ 
codecov
#TODO: replace with real black-box test
ngsmap unimplemented
./dodotest.sh
