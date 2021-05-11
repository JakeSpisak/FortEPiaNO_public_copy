from fortepianoWrapper import fortepianowrapper as fp
import numpy as np

assert fp.w_photondensity(1.2)==np.pi**2/15 * 1.2**4
print("--->Wrapper compiled and imported correctly!")
