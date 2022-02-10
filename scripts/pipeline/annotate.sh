# Draft annotation
python scripts/process/annotate.py -i data/prc/sc/integrated.h5ad -r 0.5 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -o figures

# Remove neurons
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 7,8,10,11,17,19 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 1 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 1 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 0.5 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 1 -o figures

# Remove sample specific cluster and Atros from Grey matter
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 10,16 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 0.5 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 2 -o figures

# Remove possible Astro-Oligo doublets
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 14 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 0.5 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 3 -o figures