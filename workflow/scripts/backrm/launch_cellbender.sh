module load system/singularity/3.9.2

#######################
# Important variables #
path_to_samples="data/raw/sn/"
model=full
expected_cells=8000
total_droplets_included=50000
fpr=0.01
epochs=150
posterior_batch_size=5
cells_posterior_reg_calc=50
#                     #
#######################
mkdir -p logs/run_cellbender/
for path_to_sample in ${path_to_samples}*; do
    sample_id=$(basename $path_to_sample)
    echo $sample_id
    mkdir -p data/prc/sn/$sample_id/
    inp_path=data/raw/sn/$sample_id/raw_feature_bc_matrix.h5
    out_path=data/prc/sn/$sample_id/cellbender_matrix.h5
    sbatch --output=logs/run_cellbender/$sample_id \
    --export=inp_path=$inp_path,out_path=$out_path,model=$model,expected_cells=$expected_cells,total_droplets_included=$total_droplets_included,fpr=$fpr,epochs=$epochs,posterior_batch_size=$posterior_batch_size,cells_posterior_reg_calc=$cells_posterior_reg_calc workflow/scripts/backrm/sbatch_cellbender.sh
done
