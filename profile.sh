#!/bin/bash
gprof bin/nuDens.exe | gprof2dot | dot -Tpng -o output.png
