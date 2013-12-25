#!/bin/bash

./psf_2D.py CalMomPass8 -b
./psf_2D.py CalMomPass7 -b
mv *.pdf figs/.
mv *.png figs/.
