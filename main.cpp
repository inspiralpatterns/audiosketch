#include <iostream>
#include <cmath>
#include "WavFile.h"

using namespace std;

#define PI 3.14159265
int main() {
    cout << "Variable time delay with circular buffer" << endl;
    int dtime = 50;                                // delay time is in ms
    int var_dt = 8;                                // variable delay time in ms

    // in file
    string filename("files/Vox.wav");
    WavInFile input(filename.c_str());
    int sr = input.getSampleRate();
    int bits = input.getNumBits();
    int nchan = input.getNumChannels();
    int samples = input.getNumSamples();

    float dt = dtime*float(sr)/1000;            // delay in samples
    float vdt = var_dt*float(sr)/1000;          // variable delay in samples

    float gain = 0.0;                           // feedback gain
    float wet = 0.4;                            // out balance
    float dry = 1 - wet;

    // out file
    WavOutFile out("files/out.wav", sr, bits, nchan);

    // memory allocation for input and output
    auto *in_data = new float[samples+int(dt)];
    auto *out_data = new float[samples+int(dt)];
    cout << samples << ' ' << bits << ' ' << nchan << endl;
    auto *del = new float[int(dt)];
    std::fill(in_data, in_data+samples+int(dt), 0);
    std::fill(out_data, out_data+samples+int(dt), 0);
    std::fill(del, del+int(dt), 0);

    float *rptr = del;                          // reading pointer
    float *wptr = del;                          // writing pointer
    float *mptr;                                // pointer for modulation

    // modulation lookup table
    auto *table = new float[sr];
    for (int j = 0; j < sr; ++j){
        table[j] = vdt/2 * (1 - (float) cos(2 * PI * j/(sr)));
    }

    // read input file
    input.read(in_data, samples);
    float frac;                                 // fractional delay part
    float interp_sample;

    // delay section
    int j = 0;
    for (int i = 0; i < samples+int(dt); i++) {
        // read operations
        frac = table[j] - floor(table[j]);
        mptr = rptr - (int) table[j];
        if (mptr < del) {mptr += (int) dt;}
        if (mptr >= &del[(int)dt]) {mptr -= (int) dt;}
        // variable time delay ----> y(n) = x(n) + ax(n - d(n))
        interp_sample = (1-frac) * (*mptr) + frac * *(mptr-1 < del? mptr-1+(int) dt: mptr-1);
        out_data[i] = wet * interp_sample + dry * in_data[i];
        out_data[i] /= 2;
        // write
        *wptr = in_data[i] + gain * out_data[i];
        // update pointer
        ++rptr;
        ++ wptr;
        if (wptr == &del[(int)dt]) {wptr = del;}
        if (rptr == &del[(int)dt]) {rptr = del;}
        j = (j == (sr)? 0 : ++j);
    }

    // write out file
    out.write(out_data, samples+int(dt));

    delete[] in_data;
    delete[] out_data;
    delete[] del;

    return 0;
}