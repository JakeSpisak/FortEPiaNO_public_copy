#!/bin/bash
python -c 'import matplotlib;matplotlib.use("Qt5Agg");import matplotlib.pyplot as plt;from nuDensOutput import NuDensRun;run = NuDensRun("'$1'", nnu='$2', plots='$3')'
