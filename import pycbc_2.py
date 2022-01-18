import pycbc
import matplotlib.pyplot as plt
from pycbc import frame
import pylab
from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass

from gwpy.timeseries import TimeSeries
strain = TimeSeries.read("BBH1-H1.gwf","H1:HWINJ_INJECTED", start=1240642018,end=1240643042)
conditioned = strain.crop(2, 2)

qspecgram = strain.q_transform(qrange=(10,90),frange=(30, 1024),outseg=(1240642018,1240643042))
plot=qspecgram.plot(ﬁgsize=[16,8])
ax=plot.gca()
ax.set_xscale('seconds')
ax.set_yscale('log')
ax.set_ylim(30,1024)
ax.set_xlim(1240642018-1.,1240643042+0.5)
ax.set_ylabel('Frequency [Hz]')
ax.grid(True,axis='y', which='both')
ax.colorbar(cmap='viridis',label='Normalized energy')
plot.saveﬁg('BBH1-H1_1.png')
