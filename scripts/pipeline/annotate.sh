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

# Remove sample specific clusters and Atros from Grey matter
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 10,15,16 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 0.5 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 2 -o figures

# Remove possible Astro-Oligo doublets
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 13,14 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 0.75 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 3 -o figures

# Remove possible Micro-Oligo doublets
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 13 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 4 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 4 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 0.75 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 4 -o figures

# Final Annotation
python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 0.75 -o data/prc/sc/ -d 0:Oligos,1:Astros,2:Oligos,3:Oligos,4:Oligos,5:Microglia,6:Oligos,7:OPC,8:Endothelia,9:T-cells,10:Astros_r,11:Astros_c,12:B-cells,13:Macrophages_f,14:Stroma,15:Oligos_d,16:Astros
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 5 -o figures