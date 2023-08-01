#!/bin/bash

# bash loop (for all combinations of interest)
for output in cellbender; do
   for deconv_model in all condition; do
       for dc_model in ulm wmean; do   
          for ct_metric in abunds props_ilr; do
              for image_features in None histogram summary texture; do
                  for n_factors in 10; do
                      python scripts/process/run_mofa.py --output $output --deconv_model $deconv_model --dc_model $dc_model --ct_metric $ct_metric --image_features $image_features --n_factors $n_factors --recompute True
                  done
              done
          done
      done
   done
done