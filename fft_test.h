
/*  fft_test.h
    Test for external spectral processing library
    Original code:  Carmine-Emanuele Cella
    Edit:           Mattia Paterna
*/

#ifndef FFT_TEST_H
#define FFT_TEST_H

#include <stdexcept>
#include <stdlib.h>
#include <iostream>
#include <cmath>

using namespace std;

const double TWOPI = M_PI * 2.;

// dsp parameters
const int fftlen = 1024;      // fft length
const int olap = 512;         // overlap-add
const int outhop = fftlen / olap;
const int hopsize = (int) ((float) outhop / 1.0);
const float scale = 1;        // output scaling factor
const float tstretch = 1.0;
// memory allocation
static float* in_buffer = new float[fftlen];
static float* out_buffer = new float[fftlen];
static float* buffer = new float[fftlen];
static float* amp = new float[fftlen];
static float* freq = new float[fftlen];
static float* ampsyn = new float[fftlen];
static float* freqsyn = new float[fftlen];
static float* old_phi = new float[fftlen];
static float* old_phisyn = new float[fftlen];
static float* workspace = new float[fftlen * 2];
static int outsamples = (int) (fftlen * 1.0); // safety
static float* output_data = new float[outsamples];
// guard point to wrap pointer to workspace
static float* workspace_guard_point = workspace;
static float* out_guard_point = output_data;

static int j = 0;

// make envelope window, according to given coefficients
float * makeWindow (int N, float a0, float a1, float a2) {
  float* out = new float[fftlen];
    // .5, .5, 0     --> hanning
    // .54, .46, 0   --> hamming
    // .42, .5, 0.08 --> blackman
    for (int i = 0; i < N; ++i) {
        out[i] = a0 - a1 * cos ((TWOPI * (float) i) / (N - 1)) + a2 *
        	cos ((2 * TWOPI * (float) i) / (N - 1)); // hamming, hann or blackman
    }

  return out;
}

// make envelope window
static float * window = makeWindow(fftlen, 0.5, 0.5, 0.0);
static float* window_guard_point = window;


// compute fft and ifft(as fft forward with inverse sign)
void fft (float *fftBuffer, long fftFrameSize, long sign) {
    float wr, wi, arg, *p1, *p2, temp;
    float tr, ti, ur, ui, *p1r, *p1i, *p2r, *p2i;
    long i, bitm, j, le, le2, k;

    for (i = 2; i < 2*fftFrameSize-2; i += 2) {
        for (bitm = 2, j = 0; bitm < 2*fftFrameSize; bitm <<= 1) {
            if (i & bitm) j++;
            j <<= 1;
        }
        if (i < j) {
            p1 = fftBuffer+i;
            p2 = fftBuffer+j;
            temp = *p1;
            *(p1++) = *p2;
            *(p2++) = temp;
            temp = *p1;
            *p1 = *p2;
            *p2 = temp;
        }
    }
    for (k = 0, le = 2; k < (long)(log(fftFrameSize)/log(2.)); k++) {
        le <<= 1;
        le2 = le>>1;
        ur = 1.0;
        ui = 0.0;
        arg = M_PI / (le2>>1);
        wr = cos(arg);
        wi = sign*sin(arg);
        for (j = 0; j < le2; j += 2) {
            p1r = fftBuffer+j;
            p1i = p1r+1;
            p2r = p1r+le2;
            p2i = p2r+1;
            for (i = j; i < 2*fftFrameSize; i += le) {
                tr = *p2r * ur - *p2i * ui;
                ti = *p2r * ui + *p2i * ur;
                *p2r = *p1r - tr;
                *p2i = *p1i - ti;
                *p1r += tr;
                *p1i += ti;
                p1r += le;
                p1i += le;
                p2r += le;
                p2i += le;
            }
            tr = ur*wr - ui*wi;
            ui = ur*wi + ui*wr;
            ur = tr;
        }
    }
}

// rectangular to polar coordinates conversion (required after fft)
void convert (const float* cbuffer, float* amp, float* freq, float* oldPhi,
              int N, int hop, float R) {
    float osamp = N / hop;
    float freqPerBin = R / (float) N;
    float expct = TWOPI * (float) hop / (float) N;
    for (int i = 0; i < N; ++i) {
        amp[i] = 2. * sqrt (cbuffer[2 * i] * cbuffer[2 * i]
                            + cbuffer[2 * i + 1] * cbuffer[2 * i + 1]);
        float phase = atan2 (cbuffer[2 * i + 1], cbuffer[2 * i]);
        float tmp = phase - oldPhi[i];
        oldPhi[i] = phase;
        tmp -= (float) i * expct;
        int qpd = (int) (tmp / M_PI);
        if (qpd >= 0) qpd += qpd & 1;
        else qpd -= qpd & 1;
        tmp -= M_PI * (float) qpd;
        tmp = osamp * tmp / (2. * M_PI);
        tmp = (float) i * freqPerBin + tmp * freqPerBin;
        freq[i] = tmp;
    }
}

// polar to rectangular coordinates conversion (required before inverse fft)
void unconvert (const float* amp, const float* freq, float* oldPhi,
                float* cbuffer, int N, int hop, float R) {
    float osamp = N / hop;
    float freqPerBin = R / (float) N;
    float expct = TWOPI * (float) hop / (float) N;

    for (int i = 0; i < N; ++i) {
        float tmp = freq[i];

        tmp -= (float)i * freqPerBin;
        tmp /= freqPerBin;
        tmp = TWOPI * tmp / osamp;
        tmp += (float) i * expct;
        oldPhi[i] += tmp;

        float phase = (oldPhi[i]); // - oldPhi[i - 1] - oldPhi[i + 1]);

        cbuffer[2 * i] = amp[i] * cos (phase);
        cbuffer[2 * i + 1] = amp[i] * sin (phase);
    }
}

// compute spectral processing (including fft forward and inverse)
void process (float* workspace,
	float* amp, float* freq, float* old_phi, float* old_phisyn,
	float* ampsyn, float* freqsyn, float pitch,
  int fftlen, int hop, int outhop, int sr) {
	  fft (workspace, fftlen, -1); // forward fft
    convert (workspace, amp, freq, old_phi, fftlen, hop, sr);
    // anti-noise - obs: fixed treshold
    float threshold = 0.1;
    for (int i = 0; i < fftlen; ++i) {
        if (amp[i] < threshold) amp[i] = 0;
    }
    // processing spectra
    memset (ampsyn, 0, sizeof (float) * fftlen);
    memset (freqsyn, 0, sizeof (float) * fftlen);
    int hspect = fftlen * pitch;
    for (int i = 0; i < fftlen; ++i) {
        int idx = ((int) i / pitch);
        if (idx < hspect) {
            ampsyn[i] += amp[idx];
            freqsyn[i] = freq[idx] * pitch;
        }
    }
  unconvert (ampsyn, freqsyn, old_phisyn, workspace, fftlen, outhop, sr);
	fft (workspace, fftlen, 1); // inverse fft
}

/*  to be executed at sample-rate in Faust
    input argument from Faust:
    - in_sample: single input sample
    - sr
    - ps (pitch shifting ratio)
    output argument:
    - out_sample: single output sample
    obs: tstretch fixed for now.
*/
float dsp_process(float in_sample, int sr, float ps){
  // heuristic for normalization
  float norm = (1. / (fftlen * olap));
  if (tstretch != 1) norm *= (tstretch * 1.2);
  if (ps > 1) norm = pow (norm, 1. / (ps * .6));
  else if (ps < 1) norm = pow (norm, ps * 1.75);
  norm *= scale;
  /*
  // version without overlap-add
  *workspace++ = in_sample * (*window++);
  *workspace++ = 0;

  //if pointer at end of workspace
  if (workspace == workspace_guard_point+(fftlen*2)) {
    // wrap workspace and window pointers
    workspace = workspace_guard_point;
    process (workspace, amp, freq, old_phi, old_phisyn, ampsyn,
            freqsyn, ps, fftlen, hopsize, outhop, sr);
    // write on output buffer
    output_data = out_guard_point;
    window = window_guard_point;
    for (int k = 0; k < fftlen; k++) {
      output_data[k] = (workspace[2 * k] * norm * window[k]);
    }
  }
  float out_sample = *output_data++;
  */

  // version with overlap-add
  j %= fftlen;
  // compute FFT every hopsize
  if (j % olap == 0){
    // fill workspace
    for (int k = j; k < j + fftlen; ++k) {
      workspace[2 * (k - j)] = in_buffer[k % fftlen] * window[k - j];
      workspace[2 * (k - j) + 1] = 0;
    }
    // spectral processing
    process (workspace, amp, freq, old_phi, old_phisyn, ampsyn,
      freqsyn, ps, fftlen, hopsize, outhop, sr);
    // overlap-add
   for (int k = j; k < j + fftlen; ++k) {
      out_buffer[k % fftlen] = (workspace[2 * (k - j)]
      * window[k - j] * norm);
    }
  }

  in_buffer[j] = in_sample * window[j];
  float out_sample = out_buffer[j];
  cout << out_sample << endl;
  j += 1;

  return out_sample;
}

#endif
