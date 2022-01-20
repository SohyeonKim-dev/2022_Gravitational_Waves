import pycbc
import matplotlib.pyplot as plt
from pycbc import frame
import pylab
from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass
from gwpy.timeseries import TimeSeries
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.waveform import get_td_waveform
from pycbc.filter import matched_filter
from pycbc.vetoes import power_chisq
from pycbc.events.ranking import newsnr 
import numpy

strain = frame.read_frame("BBH2-H1.gwf","H1:HWINJ_INJECTED")

strain = highpass(strain, 30.0)
#strain = resample_to_delta_t(strain, 1.0/2048)
conditioned = strain.crop(2, 2)

psd = conditioned.psd(4)
psd = interpolate(psd, conditioned.delta_f)
psd = inverse_spectrum_truncation(psd, int(4 * conditioned.sample_rate),low_frequency_cutoff=30)

for i in range(90):
    m1=i+10
    for j in range(90):
        m2=j+10
        if m1>m2:
            print(m1)
            hp, hc = get_td_waveform(approximant="SEOBNRv4_opt", mass1=m1, mass2=m2, delta_t=conditioned.delta_t, f_lower=30)
            hp.resize(len(conditioned))
            template = hp.cyclic_time_shift(hp.start_time)

            snr = matched_filter(template, conditioned, psd=psd, low_frequency_cutoff=30)
            snr = snr.crop(4 + 4, 4)

            nbins=26
            chisq = power_chisq(hp, conditioned, nbins, psd, low_frequency_cutoff=30.0)
            chisq = chisq.crop(4+4,4)
            dof = nbins * 2 - 2
            chisq /= dof
            nsnr = newsnr(abs(snr), chisq)
            peak = nsnr.argmax()
            snrp = nsnr[peak]
            time = snr.sample_times[peak]
            if abs(snrp) > 10:
                print("We found a signal at {}s with SNR {} m1:{} m2:{}".format(time, abs(snrp), m1, m2))
