#include <iostream>
#include <cmath>
#include <random>
#include "WavFile.h"

using namespace std;

#define PI 3.14159265
int main() {
    cout << "Delay line with circular buffer" << endl;
    int dtime = 500;                               // delay time is in ms

    // in file
    string filename("files/cello.wav");
    WavInFile input(filename.c_str());
    int sr = input.getSampleRate();
    int bits = input.getNumBits();
    int nchan = input.getNumChannels();
    int samples = input.getNumSamples();

    float dt = dtime*float(sr)/1000;            // delay in samples
    float gain = 0.3;                           // feedback gain
    float wet = 0.7;                            // out balance
    float dry = 1 - wet;

    // out file
    WavOutFile out("files/out.wav", sr, bits, nchan);

    // memory allocation for input and output
    auto *in_data = new double[samples+int(dt)];
    auto *out_data = new double[samples+int(dt)];
    auto *del = new double[int(dt)];
    std::fill(in_data, in_data+samples+int(dt), 0);
    std::fill(out_data, out_data+samples+int(dt), 0);
    std::fill(del, del+int(dt), 0);

    double *rptr, *wptr = del;

    // LFO wavetable: fill with hamming window
    auto *lfo_window = new double[sr];
    for (int j = 0; j < sr; ++j) {
        lfo_window[j] = 0.5 * (1 - cos(2 * PI * j/(sr)));
    }

    // read input file
    input.read(in_data, samples);

    // initialise generator to sample for jittering part
    default_random_engine generator;
    normal_distribution<float> distribution(0, 1);
    double sample, x_jit = 0;
    double y_jit = 50;
    double lfo_val, x_in;

    // delay section
    int j = 0;
    for (int i = 0; i < samples+int(dt); i++) {
        // sample from normal distribution (100 Hz)
        // todo: interpolation missing
        if (i % (int)floor(sr/10) == 0) {
            sample = distribution(generator);
            if ((sample >= -1) && (sample <= 1)) {
                // keep the value only
                x_jit = sample;
            }
        }

        // apply jittering to LFO wavetable
        lfo_val = lfo_window[j] + x_jit;
        lfo_val *= (y_jit * 0.1);
        lfo_val += 1 - (y_jit*0.1);
        lfo_val *= lfo_window[j];

        // apply LFO to in source
        x_in = in_data[i] * lfo_val;

        // read operations
        rptr = wptr - (int) dt;
        if (rptr < del) {rptr += (int) dt;}
        if (rptr >= &del[(int)dt]) {rptr -= (int) dt;}
        out_data[i] = *rptr;

        // sum original and fx
        out_data[i] *= wet;
        out_data[i] += dry * in_data[i];
        // hard clipping
        out_data[i] = out_data[i] < -1 ? -1 : out_data[i];
        out_data[i] = out_data[i] > 1 ? 1 : out_data[i];

        // write
        *wptr = x_in + gain * out_data[i];

        // update pointer
        ++rptr;
        ++ wptr;
        ++j;
        if (wptr == &del[(int)dt]) {wptr = del;}
        if (rptr == &del[(int)dt]) {rptr = del;}
        if (j == sr-1) {j = 0;}
    }

    // write out file
    out.write(out_data, samples+int(dt));

    delete[] in_data;
    delete[] out_data;
    delete[] del;
    delete[] lfo_window;

    return 0;
}