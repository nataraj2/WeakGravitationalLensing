rm -rf out
mpicxx -I/usr/local/include -L/usr/local/lib -lfftw3 power_spectrum.cpp -o out
./out
