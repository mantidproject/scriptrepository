"""
This is a script that  produces a 2x2 figure containing sinusoidal waveforms, 
two of which have error bars. this is another bit of doc rext to see if it updates
testing autoupdate
Another change so it belongs to me
"""


# The following line helps with future compatibility with Python 3
# print must now be used as a function, e.g print('Hello','World')
from __future__ import (absolute_import, division, print_function, unicode_literals)

# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *

import matplotlib.pyplot as plt

import numpy as np
import scipy.stats as stats
import math

plt.close('all')

fig = plt.figure(figsize=[9.0, 4.7766])
for i in range(4):
    ax = fig.add_subplot(2, 2, i+1, projection='mantid')
    ax.set_title("Title {}".format(i))
    x = np.linspace(0, 10, 100)
    y = (i+1)*np.sin(x + i)
    e = np.sqrt(np.abs(y))
    ws = CreateWorkspace(x, y, DataE=e, OutputWorkspace="ws"+str(i))
    if i % 2 == 1:
        ax.errorbar(ws, label=str(i), errorevery=2, elinewidth=0.6, specNum=1)
    else:
        ax.plot(ws, label=str(i), specNum=1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    ax.legend()
plt.show()