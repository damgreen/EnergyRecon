#!/bin/bash

./psf_2D.py CalMomP8 -b
./psf_2D.py CalMomP7 -b
mv *.pdf figs/.
mv *.png figs/.
