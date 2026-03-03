#!/bin/bash
# run_2025pht_models.sh

# Every linear combination of models to run in case they need to be regenerated

python ../analysis.py 2025pht rsg --dust-comp grp --nsteps 25000 --use-variance --dust-width 2
python ../analysis.py 2025pht rsg --dust-comp sil --nsteps 25000 --use-variance --dust-width 2
python ../analysis.py 2025pht-early rsg --dust-comp grp --nsteps 25000 --use-variance --dust-width 2
python ../analysis.py 2025pht-early rsg --dust-comp sil --nsteps 25000 --use-variance --dust-width 2
python ../analysis.py 2025pht-short rsg --dust-comp grp --nsteps 25000 --use-variance --dust-width 2
python ../analysis.py 2025pht-short rsg --dust-comp sil --nsteps 25000 --use-variance --dust-width 2

python ../analysis.py 2025pht rsg --dust-comp grp --nsteps 25000 --use-variance --dust-width 10
python ../analysis.py 2025pht rsg --dust-comp sil --nsteps 25000 --use-variance --dust-width 10
python ../analysis.py 2025pht-early rsg --dust-comp grp --nsteps 25000 --use-variance --dust-width 10
python ../analysis.py 2025pht-early rsg --dust-comp sil --nsteps 25000 --use-variance --dust-width 10
python ../analysis.py 2025pht-short rsg --dust-comp grp --nsteps 25000 --use-variance --dust-width 10
python ../analysis.py 2025pht-short rsg --dust-comp sil --nsteps 25000 --use-variance --dust-width 10

python ../analysis.py 2025pht blackbody --nsteps 5000 --use-variance

python ../plotting.py 2025pht sed rsg,blackbody --dust-comp grp,sil --use-variance --dust-width 2
python ../plotting.py 2025pht-early sed rsg --dust-comp grp,sil --use-variance --dust-width 2 --plot-title "Only 5 Feb 2024"
python ../plotting.py 2025pht-short sed rsg --dust-comp grp,sil --use-variance --dust-width 2 --plot-title "Only NIRCam"
python ../plotting.py 2025pht sed rsg --dust-comp grp,sil --use-variance --dust-width 10 --plot-title r"$R_{\rm out}/R_{\rm in}$=10"

