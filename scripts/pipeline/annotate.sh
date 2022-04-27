echo '## Annotating cell types ##'
# Draft annotation
python scripts/process/annotate.py -i data/prc/sc/integrated.h5ad -r 0.5 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -o figures

# Remove neurons and doublet clusters
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 7,9,10,12,13,17 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 1 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 1 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.25 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 1 -o figures

# Remove stressed cells and Astros from Grey Matter
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 17,25 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 2 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.25 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 2 -o figures

# Remove cluster of doublets
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 21 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 3 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 2.0 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 3 -o figures

# Remove cluster of microglia/oligos doublets
python scripts/process/remerge.py -i data/prc/sc/annotated.h5ad -c 28 -o data/prc/sc/
python scripts/plot/merge_metrics.py -i data/prc/sc/annotated.h5ad -v 4 -o figures/

python scripts/process/reintegrate.py -i data/prc/sc/annotated.h5ad -o data/prc/sc/
python scripts/plot/integrate_metrics.py -i data/prc/sc/annotated.h5ad -v 4 -o figures/

python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.0 -o data/prc/sc/
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 4 -o figures

# Final Annotation
python scripts/process/annotate.py -i data/prc/sc/annotated.h5ad -r 1.0 -o data/prc/sc/ -d 0:Oligos,1:Oligos,2:Oligos,3:Microglia,4:Astros,5:Oligos,6:Astros,7:Oligos,8:OPC,9:Endothelia,10:Oligos,11:Oligos,12:T-cells,13:Astros,14:Astros_r,15:Astros_c,16:B-cells,17:Macrophages_f,18:Stroma,19:Astros,20:Oligos_d
python scripts/plot/annotate_metrics.py -i data/prc/sc/annotated.h5ad -m data/markers.csv -v 5 -o figures
