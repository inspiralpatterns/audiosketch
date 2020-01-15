import("stdfaust.lib");
import("analyzers.lib");


sine = os.oscsin(440.0);
N = 1024;
process = sine : an.rtocv(N) : an.fft(N);
