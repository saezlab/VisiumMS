#######################
# Important variables #
path_to_regres="/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/data/prc/visium/reg_model/inf_aver.csv"
path_to_slides="/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/data/raw/visium/"
path_to_output="/net/data.isilon/ag-saez/bq_pbadia/VisiumMS/data/prc/visium/"
#                     #
#######################


for path_to_slide in ${path_to_slides}*; do
    sbatch --export=path_to_regres=$path_to_regres,path_to_slide=$path_to_slide,path_to_output=$path_to_output \
    -J dcnv_$(basename $path_to_slide) -o logs/dcnv_$(basename $path_to_slide).out scripts/launchers/deconv_slide.sh
done
