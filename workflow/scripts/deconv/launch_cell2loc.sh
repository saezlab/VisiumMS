module load system/singularity/3.9.2

#######################
# Important variables #
path_to_samples='data/prc/vs/'
reg_path='data/prc/vm_regmodel.csv'
n_cells_spot=5
d_alpha=20
max_epochs=30000
res_path='results/deconv/'
#                     #
#######################
mkdir -p logs/run_cell2loc/
mkdir -p $res_path
for path_to_sample in ${path_to_samples}*; do
    sample_id=$(basename $path_to_sample)
    echo $sample_id
    slide_path=$path_to_samples/$sample_id/adata.h5ad
    plot_path=$res_path/train_$sample_id.pdf
    out_path=$path_to_samples/$sample_id
    sbatch --output=logs/run_cell2loc/$sample_id \
    --export=slide_path=$slide_path,reg_path=$reg_path,n_cells_spot=$n_cells_spot,d_alpha=$d_alpha,max_epochs=$max_epochs,plot_path=$plot_path,out_path=$out_path workflow/scripts/deconv/sbatch_deconv.sh
done
