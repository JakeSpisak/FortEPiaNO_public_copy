import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import numpy as np
from nuDensOutput import colors, styles, stripRepeated, NuDensRun

fNO = NuDensRun("OUT/3nu_20f_xc_no/", label="NO")
fIO = NuDensRun("OUT/3nu_20f_xc_io/", label="IO")
fNOd = NuDensRun("OUT/3nu_20f_d_no/", label="NO damp")
fNOd = NuDensRun("OUT/3nu_20f_d_io/", label="IO damp")
fnoosc = NuDensRun("OUT/3nu_20f_od_no/", label="no osc", full=False)
