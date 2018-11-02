import numpy as np
import matplotlib.pyplot as plt

def plot_tp(filt, cutoff, label, name):
     wave = filt['wave']*1e4
     pce = filt['pce']*1.5
     wave_use = wave[pce > cutoff]
     pce_use = pce[pce > cutoff]
     w_min = np.min(wave_use)
     w_max = np.max(wave_use)
     med = np.median(pce_use)
     mean = np.mean(pce_use)
     mx = np.max(pce_use)
     plt.clf()
     plt.plot(wave, pce, linestyle='--', color=PANK)
     plt.plot(wave_use, pce_use, color=DBLUE)
     plt.axhline(mx, linestyle=':', color=LPURPLE, label='max')
     plt.axhline(med, linestyle=':', color=LBLUE, label='med')
     plt.axhline(mean, linestyle=':', color=DGREEN, label='mean')
     plt.legend(bbox_to_anchor=(.3,.25), loc=2, borderaxespad=0.)
     plt.ylabel(label + ' throughput')
     plt.xlabel('wavelength [angstrom]')
     plt.savefig(name + '.png')
     print("filters['" + label + "'] = {'min_max': (" + str(w_min) + ", " + str(w_max) + "), 'max': " +str(mx) +  ", 'med': " + str(med) + ", 'mean': " + str(mean) + "}")

