
annotation_dict = {
    "cellbender": {
        "0.75": {
            '0':'Oligos',
            '1':'Oligos',
            '2':'Astros',
            '3':'Oligos',
            '4':'Microglia',
            '5':'Oligos',
            '6':'Oligos', 
            '7':'Oligos', # lots of MT genes, could it be doublets? 
            '8':'OPC', # inside this cluster there are the Oligo differentiating but at the moment we cant separate them
            '9':'Neurons',
            '10':'Oligos',
            '11':'Endothelia',
            '12':'Neurons', 
            '13':'Neurons',
            '14':'Neurons',
            '15':'T_cells',
            '16':'Astros', # this is a new cluster of astros that we did not have before and its in both datasets INTERESTING
            '17':'B_cells',
            '18':'Astros_c',
            '19':'Stroma',
            '20':'Neurons',
            '21':'Random_cluster' #This is the cluster that i mentioned that i would delete
        }
    },
    "cellranger": {
        "0.75": { 
        '0':'Oligos',
        '1':'Oligos',
        '2':'Oligos',
        '3':'Microglia',
        '4':'Astros',
        '5':'DOUBLETS',  #REMOVE: mix oligo, astro, b cell,endothelia
        '6':'Oligos',
        '7':'OPC', #inside this cluster there are as well differentiating oligos
        '8':'Astros_r', 
        '9':'Neurons',
        '10':'Neurons',
        '11':'Endothelia',
        '12':'Neurons',
        '13':'T_cells',
        '14':'Neurons', 
        '15':'Macrophages', # Super inflammed and a lot of ribosomal genes
        '16':'Astros', # the new astros cluster
        '17':'B_cells',
        '18':'Astros_weird', # It has some astros markers, but also a lot of undefined genes. Would most probably REMVOVE
        '19':'Astros_c',
        '20':'Stroma',
        '21':'Neurons',
        '22':'Neurons',
        '23':'Oligos',
        '24':'Random_cluster' # REMOVE
        }
    }
}
