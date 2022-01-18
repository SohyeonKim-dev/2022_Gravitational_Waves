import pycbc
import matplotlib.pyplot as plt
from pycbc import frame
import pylab
from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass
from gwpy.timeseries import TimeSeries

strain = frame.read_frame("BBH1-H1.gwf","H1:HWINJ_INJECTED")

strain = highpass(strain, 30.0)
#strain = resample_to_delta_t(strain, 1.0/2048)
conditioned = strain.crop(2, 2)

from pycbc.psd import interpolate, inverse_spectrum_truncation
psd = conditioned.psd(4)
psd = interpolate(psd, conditioned.delta_f)
psd = inverse_spectrum_truncation(psd, int(4 * conditioned.sample_rate),low_frequency_cutoff=30)

from pycbc.waveform import get_td_waveform
m = 36 
hp, hc = get_td_waveform(approximant="SEOBNRv4_opt", mass1=m, mass2=m, delta_t=conditioned.delta_t, f_lower=30)
hp.resize(len(conditioned))
pylab.figure()
pylab.title('Before shifting')
pylab.plot(hp.sample_times, hp)
pylab.xlabel('Time (s)')
pylab.ylabel('Strain')
pylab.savefig('1st.png')
template = hp.cyclic_time_shift(hp.start_time)
pylab.figure()
pylab.title('After shifting')
pylab.plot(template.sample_times, template)
pylab.xlabel('Time (s)')
pylab.ylabel('Strain')
pylab.savefig('2nd.png')

from pycbc.filter import matched_filter
import numpy
snr = matched_filter(template, conditioned, psd=psd, low_frequency_cutoff=30)
snr = snr.crop(4 + 4, 4)
pylab.figure(figsize=[10, 4])
pylab.plot(snr.sample_times, abs(snr))
pylab.ylabel('Signal-to-noise')
pylab.xlabel('Time (s)')
pylab.savefig("victory.png")

peak = abs(snr).numpy().argmax()
snrp = snr[peak]
time = snr.sample_times[peak]

print("We found a signal at {}s with SNR {}".format(time, abs(snrp)))

from pycbc.filter import sigma
dt = time - conditioned.start_time
aligned = template.cyclic_time_shift(dt)

aligned /= sigma(aligned, psd=psd, low_frequency_cutoff=30.0)

aligned = (aligned.to_frequencyseries() * snrp).to_timeseries()
aligned.start_time = conditioned.start_time

white_data = (conditioned.to_frequencyseries() / psd**0.5).to_timeseries()
white_template = (aligned.to_frequencyseries() / psd**0.5).to_timeseries()

white_data = white_data.highpass_fir(30., 512).lowpass_fir(300, 512)
white_template = white_template.highpass_fir(30, 512).lowpass_fir(300, 512)

white_data = white_data.time_slice(1240642818.0134277-.2, 1240642818.0134277+.1)
white_template = white_template.time_slice(1240642818.0134277-.2, 1240642818.0134277+.1)

pylab.figure(figsize=[15, 3])
pylab.plot(white_data.sample_times, white_data, label="Data")
pylab.plot(white_template.sample_times, white_template, label="Template")
pylab.legend()
pylab.savefig('victory2.png')

subtracted = conditioned - aligned

for data, title in [(conditioned, 'Original H1 Data'),
                    (subtracted, 'Signal Subtracted from H1 Data')]:

    t, f, p = data.whiten(4, 4).qtransform(.001, logfsteps=100, qrange=(8, 8), frange=(20, 512))
    pylab.figure(figsize=[15, 3])
    pylab.title(title)
    pylab.pcolormesh(t, f, p**0.5, vmin=1, vmax=6)
    pylab.yscale('log')
    pylab.xlabel('Time (s)')
    pylab.ylabel('Frequency (Hz)')
    pylab.xlim(1240642818.0134277 - 2, 1240642818.0134277 + 1)
    pylab.show()
    pylab.savefig("victory3.png")
