#!/bin/bash

# ����·��
cd /home/wzy/QES-Public-main/buildGPU/

./qesPlume/qesPlume  -q /home/wzy/QES-Public-main/LMZ_QES_CESTA1987/CESTA1987plume.xml   -w /home/wzy/QES-Public-main/LMZ_QES_CESTA1987/Output/lmz_windsWk.nc   -t  /home/wzy/QES-Public-main/LMZ_QES_CESTA1987/Output/lmz_turbOut.nc   -o /home/wzy/QES-Public-main/LMZ_QES_CESTA1987/Output/lmz
