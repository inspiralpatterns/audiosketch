#include <iostream>
#include <cmath>
#include <random>
#include "WavFile.h"

using namespace std;

#define PI 3.14159265


/* delay function
 * input arguments:
 * wptr:        writing pointer
 * del:         delay buffer
 * dt:          delay time (samples)
 * x_in:        pointer to input buffer
 * y_out:       pointer to output buffer
 * wet/dry:     original/processed linear amplitude
 *
 * output arguments:
 * wptr: updated writing pointer to delay line
 */
double *delay(double *wptr, double *del, double dt, double x_in, double *y_out, float feedback_gain){
    // read
    double *rptr = wptr - (int) dt;
    if (rptr < del) {rptr += (int) dt;}
    if (rptr >= &del[(int)dt]) {rptr -= (int) dt;}
    *y_out = *rptr;

    // write
    *wptr = x_in + feedback_gain * (*y_out);

    // update pointer
    ++wptr;
    if (wptr == &del[(int) dt]) { wptr = del; }

    return wptr;
}


int main() {
    cout << "Oscillating | wave-based delay LFO [v 0.1]" << endl;
    int dtime = 500;    // delay time is in ms

    // in file
    string filename("files/Vox.wav");
    WavInFile input(filename.c_str());
    int sr = input.getSampleRate();
    int bits = input.getNumBits();
    int nchan = input.getNumChannels();
    int samples = input.getNumSamples();

    float dt = dtime*float(sr)/1000;    // delay in samples
    float feedback_gain = 0.1;
    float wet = 0.5;
    float dry = 1 - wet;

    // out file
    WavOutFile out("files/out.wav", sr, bits, nchan);

    // memory allocation for input and output
    auto *in_data = new double[samples+int(dt)];
    auto *out_data = new double[samples+int(dt)];
    double *origin = out_data;
    double *guard_point = &out_data[samples+int(dt)];
    auto *del = new double[int(dt)];
    std::fill(in_data, in_data+samples+int(dt), 0);
    std::fill(out_data, out_data+samples+int(dt), 0);
    std::fill(del, del+int(dt), 0);

    // LFO wavetable: fill with hamming window
    auto *wave_window = new double[sr];
    for (int j = 0; j < sr; ++j) {
        wave_window[j] = 0.5 * (1 - cos(2 * PI * j/(sr)));
    }

    // initialise generator to sample for jittering part
    default_random_engine generator;
    normal_distribution<float> distribution(0, 1);
    double x_jit;   // controls sampling freq (range 0-300)
    double y_jit;   // controls distortion (range 0-100)
    int samp_freq;
    double lfo_val, x_in;
    double a, b;    // internal jittering variables

    // read input file
    input.read(in_data, samples);

    // initialise jittering variables
    samp_freq = 2;
    x_jit = distribution(generator);
    y_jit = 100;

    // initialise pointers
    double *wptr = del;
    int j = 0;
    // iterate over input
    for (int i = 0; i < samples+int(dt); i++) {
        // sample from normal distribution - only if sampling frequency is positive
        // todo: needs interpolation
        if (samp_freq > 0) {
            if (i % (int) floor(sr / samp_freq) == 0){
                do {
                    x_jit = distribution(generator);
                } while ((x_jit < -1) || (x_jit > 1));
            }
        }

        // apply jittering to LFO wavetable
        a = wave_window[j] + x_jit;
        a *= (y_jit * 0.1);
        b = wave_window[j] * 1 - (y_jit * 0.1);
        lfo_val = (a + b) * wave_window[j];
        // apply LFO to in source
        x_in = lfo_val * in_data[i];

        // call to delay
        wptr = delay(wptr, del, dt, x_in, out_data, feedback_gain);
        // wet-dry balance
        (*out_data) *= wet;
        (*out_data) += dry * in_data[i];
        // update output buffer pointer
        ++out_data;
        if (out_data == guard_point) { cout << "guard point reached" << endl; }
        if (j == sr - 1) { j = 0; } else { ++j; }
    }

    // write out file - pointer back to origin!
    out_data = origin;
    out.write(out_data, samples+int(dt));
    cout << "file written" << endl;

    delete[] in_data;
    delete[] out_data;
    delete[] del;
    delete[] wave_window;

    return 0;
}