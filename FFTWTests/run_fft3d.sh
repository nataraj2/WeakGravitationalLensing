rm -rf out
mpicxx -I/opt/cray/pe/fftw/3.3.10.3/x86_milan/include -L/opt/cray/pe/fftw/3.3.10.3/x86_milan/lib -lfftw3 fft3d.cpp -o out
./out
