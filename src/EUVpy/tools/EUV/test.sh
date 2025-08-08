#!/bin/sh

# here is a test:
./fism2_process.py -start 20031028 -end 20031101 -euvfile euv.csv
mv fism200310_nWaves_037.dat fism200310_nWaves_037_daily.dat
mv fism200310_nWaves_037.png fism200310_nWaves_037_daily.png

# get flare data (this takes a while to download):
./fism2_process.py -start 20031028 -end 20031029 -euvfile euv.csv -flare
mv fism200310_nWaves_037.dat fism200310_nWaves_037_flare.dat
mv fism200310_nWaves_037.png fism200310_nWaves_037_flare.png
