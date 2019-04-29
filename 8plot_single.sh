#!/bin/bash
python -c 'import matplotlib;matplotlib.use("Qt5Agg");import matplotlib.pyplot as plt;from fortepianoOutput import FortEPiaNORun;run = FortEPiaNORun("'$1'", nnu='$2', plots='$3')'
