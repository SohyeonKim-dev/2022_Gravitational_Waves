import pycbc
import matplotlib.pyplot as plt
from pycbc import frame
import pylab
from pycbc.catalog import Merger
from pycbc.filter import resample_to_delta_t, highpass

# from gwpy.timeseries import TimeSeries
strain = frame.read_frame("BBH1-H1.gwf","H1:HWINJ_INJECTED")

# strain = highpass(strain, 30.0)
# strain = resample_to_delta_t(strain, 1.0/2048)
# conditioned = strain.crop(2, 2)

pylab.plot(strain.sample_times, strain)
pylab.xlabel('Time (s)')
pylab.savefig('tttt.png')
