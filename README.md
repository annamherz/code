# code

WIP , currently not fully useable

various scripts and code for version control

bash scripts - assorted bash scripts so I dont loose them
other scripts - other scripts and notebooks used
pipeline notebooks - notebooks used to run the pipeline from scratch
pipeline scripts - actual scripts used to run the FEP pipeline
python - python code of functions used in the scripts
tutorial - outline of functions and how they fit together in the overall scripts

install for running (must also have benchmark or feature-amber-fep BSS branch and its dependancies installed) using:
python setup.py install


## summary of functions

move to utils potentially ? static pickle ext, file ext, update options dict? colour dict ?

analysis
    _analysis : analyse (default options, file extension, pickle extension, validate options, set options, 
                analyse all repeats (calls normal or pickle, calc freenrg), save pickle, check overlap, plot graphs (matrix or dhdl), calculate  convergence, plot convergence, calculate truncated, plot across lambda, edgembar)
    
    _convert - static methods : into fwf, M_kcal, yml into exper dict, csv into exper dict, raw_dict into val dict, cinnabar file

    _dictionary : make_dict (comp_results, value_list_from_files, exp from fwf (val, diff), fwf net ana, cinnabar network edges, cinnabar node, experimental for network, exper form ligands)

    _network : get ligands from perts, get_info_network, get_info_network_dict, 
               network_graph (draw graph, draw ligand, draw_all_ligads, draw_perturbation, disconnected ligands, cycle_closures, cycle_closure_dict, cycle_vals, average_cycle_closures, add_weight, )

    _plotting : plotting_engines (match_dicts_to_df, prune_perturbations, set_colours, bar, scatter,            pert_dict_into_convergence_df)
                plotting_histogram (histogram, histogam distribution)

    _analysis_network : analysis_network( get_experimental, remove perturbations, remove ligands, change name, compute_results, successful perturbations, failed, drawing these and disconnected ligands, get_outliers, 
    remove outliers, sort ligands, compute_other_results, compute_convergence, compute consensus, draw graph, 
    compute cycle closures)

    _statistics : compute statistics wrapped around cinnabar

