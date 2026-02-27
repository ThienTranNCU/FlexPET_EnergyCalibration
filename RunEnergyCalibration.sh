#!/bin/bash

g++ FlexPET_EnergyCal_v2.cpp -o FlexPET_EnergyCal_v2 `root-config --cflags --glibs` -lSpectrum

./FlexPET_EnergyCal_v2 <InputFolder>

#306.82
#511







