#include <iostream>
#include "WavFile.h"
#include "fft_test.h"

using namespace std;

static int* X = new int[5];
static int * init = X;

// Test function for an understanding of pointers and their increment
void store(int n){
    // store value
    cout << X << endl;
    *X++ = n;
    cout << "value: " << *(X-1) << endl;
    // obs: addresses do not match
    if (X == init+5){
        cout << "wrapping pointer to " << init << endl;
        // X = &X[0]
        X -= 5;
    }
    cout << "---" << endl;
}


int main(int argc, char const *argv[]) {
  // in file
  string filename("files/monocello.wav");
  WavInFile input(filename.c_str());
  int sr = input.getSampleRate();
  int bits = input.getNumBits();
  int nchan = input.getNumChannels();
  int samples = input.getNumSamples();

  // out file
  WavOutFile out("files/out.wav", sr, bits, nchan);
  // memory allocation for input and output
  float *in_data = new float[samples * nchan];
  float *out_data = new float[samples * nchan];

  input.read(in_data, samples * nchan);
  for (int i = 0; i < samples * nchan; i++) {
    out_data[i] = dsp_process(in_data[i], sr, 1.0);
  }

  out.write(out_data, samples * nchan);

  return 0;
}
