echo '## Annotating cell types ##'
# Draft annotation
python scripts/process/annotate.py -i data/prc/sc/integrated.h5ad -r 0.5 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -o figures

# Remove neurons and doublet clusters
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 7,10,12,13,14,18,20,21 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 1 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 1 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.0 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 1 -o figures

# Remove sample specific clusters, Atros from Grey matter and stressed cells cluster
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 15,16,20 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.0 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 2 -o figures

# Remove cluster of neurons
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 12 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.0 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 3 -o figures

# Remove cluster of microglia/oligos doublets
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 14 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 4 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 4 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.0 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 4 -o figures

# Final Annotation
python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.0 -o data/prc/sc/ -d 0:Oligos,1:Oligos,2:Astros,3:Microglia,4:Oligos,5:Oligos,6:Oligos,7:Oligos,8:Oligos,9:OPC,10:Endothelia,11:Astros,12:Astros,13:T-cells,14:Astros_r,15:Macrophages_f,16:Astros_c,17:B-cells,18:Stroma,19:Oligos_d
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 5 -o figures
