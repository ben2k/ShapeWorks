#!/bin/bash

shapeworks readimage --name $DATA/1x2x2.nrrd translate -x 10.0 -y 10.0 -z 10.0 compare --name $DATA/translate1.nrrd

shapeworks readimage --name $DATA/1x2x2.nrrd translate -x -10.0 -y -10.0 -z -10.0 compare --name $DATA/translate2.nrrd
