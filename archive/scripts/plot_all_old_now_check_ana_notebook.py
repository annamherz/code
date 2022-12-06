###this script should plot all results for the different engines.
#first need to run analysis script.
#other input required is - 
#some way to exclude the experimental entirely if its not available?

#anna

#some way to choose which error is being plotted??

#import libraries
import BioSimSpace as BSS
import pandas as pd
import numpy as np
from numpy import log as ln
import os
import matplotlib.pyplot as plt
import itertools as it
import math
import random
import pickle
from sklearn.linear_model import LinearRegression #makes linear polynomial expression

# res_dir = '/home/anna/Documents/benchmark/1styr_report_results'
main_dir = os.environ["MAINDIRECTORY"]
res_dir = f"{main_dir}/outputs"
print(res_dir)

# TODO best to have this so gets engines and transformations from somewhere?
engine = ['AMBER','GROMACS','SOMD']
trans = ['ejm42~ejm31','ejm42~ejm55','ejm54~ejm42','ejm55~ejm54','2w~2z','67~60','60~63']
# trans = ['ejm42~ejm31','ejm42~ejm55','ejm54~ejm42','ejm55~ejm54']

# colour_dict = {"AMBER":"orange","SOMD":"darkturquoise","GROMACS":"orchid","experimental":"midnightblue"}
# # plot the convergence w time
# for tra in trans:
#     prot = "tyk2"

#     for leg in [ 'free','bound']:
#         plt.figure()
#         lines = []
#         for eng in engine:
#             # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/{leg}_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#             with open (f'/home/anna/Documents/benchmark/{prot}/pickles/{leg}_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#                 pmf_dict = pickle.load(handle)
#             lines += plt.plot(0,0,c=col, label=eng)
#             for repeat in pmf_dict:
#                 pmf = pmf_dict[repeat]
#                 x =[]
#                 y=[]
#                 yerr = []
#                 for p in pmf:
#                     x.append(p[0])
#                     y.append(p[1]*(1/BSS.Units.Energy.kcal_per_mol))
#                     yerr.append(p[2]*(1/BSS.Units.Energy.kcal_per_mol))
#                 plt.errorbar(x,y,yerr=yerr,color=col, ecolor='black')
#         plt.xlim(xmin=0,xmax=1)
#         plt.ylabel("Computed $\Delta$G$_{transformation}$ / kcal$\cdot$mol$^{-1}$")
#         plt.xlabel("Lambda")
#         labels = [l.get_label() for l in lines]
#         plt.legend(lines, labels)
#         plt.title(f"Convergence, {leg} for {tra}")
#         plt.savefig(f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/convergence_{leg}.png')

#     # plotting delta delta G

#     plt.figure()
#     lines = []
#     for eng in engine:
#         # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/bound_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#         with open (f'/home/anna/Documents/benchmark/{prot}/pickles/bound_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#             bound_pmf_dict = pickle.load(handle)
#         # with open (f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/free_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#         with open (f'/home/anna/Documents/benchmark/{prot}/pickles/free_pmf_{tra}_{eng}.pickle', 'rb') as handle:
#             free_pmf_dict = pickle.load(handle)
#         lines += plt.plot(0,0,c=colour_dict[eng], label=eng)
#         for repf,repb in zip(free_pmf_dict,bound_pmf_dict):
#             bound_pmf = bound_pmf_dict[repb]
#             free_pmf = free_pmf_dict[repf]
#             x = []
#             y = []
#             yerr = []
#             for pb,pf in zip(bound_pmf,free_pmf):
#                 x.append(pb[0])
#                 y.append((pb[1]*(1/BSS.Units.Energy.kcal_per_mol))-(pf[1]*(1/BSS.Units.Energy.kcal_per_mol)))
#                 yerr.append((pb[2]*(1/BSS.Units.Energy.kcal_per_mol))+(pf[2]*(1/BSS.Units.Energy.kcal_per_mol)))
#             plt.errorbar(x,y,yerr=yerr,color=colour_dict[eng])
#     plt.xlim(xmin=0,xmax=1)
#     plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
#     plt.xlabel("Lambda")
#     labels = [l.get_label() for l in lines]
#     plt.legend(lines, labels)
#     plt.title(f"Convergence for {tra}")
#     plt.savefig(f'/home/anna/Documents/benchmark/{prot}/outputs/{eng}/{tra}/convergence_deltadeltaG.png')
        
#open and calculate the experimental first off
print("calculating experimental transformation values...")

#ideally different analysis based on which data is from the experimental - also need to add try so doesnt fail in case not there
exper = pd.read_csv(f'{res_dir}/experimental.csv', sep=',', header='infer', index_col=None)
exper_calc = pd.DataFrame(index=['lig','Ki(microM)','deltaG(kcal/mol)'])

for lig in exper['lig']:
    Ki = exper.loc[exper['lig']== lig, 'Ki(microM)'].iloc[0]
    #calculate the deltaG in kcal/mol
    T = 297
    deltaG = ((8.314 * T * ln(Ki*(10**(-6)))) * (10**(-3))) * 0.239006
    # assuming for now standard deviation is 0.4 kcal/mol
    sd = 0.4
    exper_calc = exper_calc.append({'lig':lig, 'Ki(microM)': Ki, 'deltaG(kcal/mol)': deltaG, 'sd':sd}, ignore_index=True).dropna()

#create data frame for the considered transformations
trans_exp = pd.DataFrame(index=["lig_1", "lig_2", "freenrg", "error", "engine"])

for tra in trans:
    lig_1 = tra.split('~')[0]
    lig_2 = tra.split('~')[1]
    lig_1_enrg = exper_calc.loc[exper_calc['lig']== lig_1, 'deltaG(kcal/mol)'].iloc[0]
    lig_2_enrg = exper_calc.loc[exper_calc['lig']== lig_2, 'deltaG(kcal/mol)'].iloc[0]
    lig_1_err = exper_calc.loc[exper_calc['lig']== lig_1, 'sd'].iloc[0]
    lig_2_err = exper_calc.loc[exper_calc['lig']== lig_2, 'sd'].iloc[0]
    deltadeltaG = lig_2_enrg - lig_1_enrg
    prop_error = math.sqrt(math.pow(lig_1_err, 2) + math.pow(lig_2_err, 2))
    trans_exp = trans_exp.append({'lig_1':lig_1, 'lig_2':lig_2, 'freenrg': deltadeltaG, 'error': prop_error, 'engine':'experimental'}, ignore_index=True).dropna()

trans_exp.to_csv(f'{res_dir}/experimental_results.csv', sep=',', index=False)

#open all files with the values
exper_res = pd.read_csv(f'{res_dir}/experimental_results.csv', sep=',', header='infer', index_col=None)
somd_res = pd.read_csv(f'{res_dir}/final_summary_SOMD.csv', sep=',', header='infer', index_col=None)
amber_res = pd.read_csv(f'{res_dir}/final_summary_AMBER.csv', sep=',', header='infer', index_col=None)
gromacs_res = pd.read_csv(f'{res_dir}/final_summary_GROMACS.csv', sep=',', header='infer', index_col=None)

#all dictionaries make at start so not caught in a loop
#make dictionary for files and engine so can easily ref
eng_dict = {'AMBER' : amber_res, 'SOMD' : somd_res, 'GROMACS' : gromacs_res, 'experimental' : exper_res}
#create new dictionary for the values for each engine to keep them safe 
eng_dict_val = {'AMBER' : None, 'SOMD' : None, 'GROMACS' : None, 'experimental' : None}
#create new dictionary for the values for each engine so can later plot all together and to keep them safe 
eng_dict_plot_eng = {'AMBER' : None, 'SOMD' : None, 'GROMACS' : None}
eng_dict_plot_exp_for_eng = {'AMBER' : None, 'SOMD' : None, 'GROMACS' : None}

##calculating the values needed for linear regression
#first get values but also turn into array
x_data = pd.DataFrame(index=['transformation','freenrg','error'])
missing_x = list()

for tra in trans:
    df = eng_dict['experimental']
    lig_1 = tra.split('~')[0]
    lig_2 = tra.split('~')[1]

    try:
        if df.loc[df['lig_1'] == lig_1, 'lig_1'].iloc[0] == lig_1 and (df.loc[df['lig_2'] == lig_2, 'lig_2'].iloc[0]) == lig_2:
            deltadeltaG = df.loc[(df['lig_1'] == lig_1) & (df['lig_2'] == lig_2), 'freenrg'].iloc[0]
            error = df.loc[(df['lig_1'] == lig_1) & (df['lig_2'] == lig_2), 'error'].iloc[0]
            x_data = x_data.append({'transformation':tra, 'freenrg': deltadeltaG, 'error':error}, ignore_index=True).dropna(how='all')
    except:
        print(f'{tra} is not found in experimental')
        x_data = x_data.append({'transformation':tra, 'freenrg': 'NaN', 'error': 'NaN'}, ignore_index=True).dropna(how='all')
        missing_x.append(tra)

#write the x to the value dictionary
eng_dict_val['experimental'] = x_data
print(eng_dict_val['experimental'])

equations = pd.DataFrame(index=['engine','slope','slope_err','intercept','intercept_err','r2','r2_err','MUE','MUE_err'])
# TODO fix as this is not for written w units which is not good
for eng in engine:
    print(f'plotting and calculating for {eng}')

    y_data = pd.DataFrame(index=['transformation','freenrg','error'])
    missing_y = list()

    for tra in trans:
        df = eng_dict[eng]
        lig_1 = tra.split('~')[0]
        lig_2 = tra.split('~')[1]
        try:
            if df.loc[df['lig_1'] == lig_1, 'lig_1'].iloc[0] == lig_1 and (df.loc[df['lig_2'] == lig_2, 'lig_2'].iloc[0]) == lig_2:
                deltadeltaG = df.loc[(df['lig_1'] == lig_1) & (df['lig_2'] == lig_2), 'freenrg'].iloc[0]
                error = df.loc[(df['lig_1'] == lig_1) & (df['lig_2'] == lig_2), 'error'].iloc[0]
                y_data = y_data.append({'transformation':tra, 'freenrg': deltadeltaG, 'error':error}, ignore_index=True).dropna(how='all')
        except:
            print(f'{tra} is not found in {eng}')
            y_data = y_data.append({'transformation':tra, 'freenrg': 'NaN', 'error': 'NaN'}, ignore_index=True).dropna(how='all')
            missing_y.append(tra)

    #put into dictionary
    eng_dict_val[eng] = y_data
    print(eng_dict_val[eng])

    #now have all the data, can calculate the linear regression
    #first need to remove the NaN data from the y_data if there is any. Is okay to do now as the data frame is safely stored in the dictionary
   
    try:
        if len(missing_y) != 0:
            y_data = y_data.set_index('freenrg')
            y_data = y_data.drop('NaN', axis=0)
            y_data = y_data.reset_index()
            y_data = y_data.set_index('transformation')
            y_data = y_data.reset_index()
        elif len(missing_y) == 0:
            y_data = y_data
    except:
        pass

    #x needs to be reshaped bc it needs to be 2 dimensional, ie one column and as many rows as needed - does it at this stage tho
    x = (np.array(x_data.iloc[0:len(x_data),1])).reshape(-1,1)
    y = np.array(y_data.iloc[0:len(y_data), 1])

    #if the x and y are not of the same length, need to remove the missing values from the data so can model
    if len(x) > len(y): #if more x, need to remove the ones from x that are missing in y
        for tra in missing_y:
            x_data = x_data.set_index('transformation')
            x_data = x_data.drop(tra, axis=0)
            x_data = x_data.reset_index()
    
    elif len(x) < len(y): #if more y, need to remove the ones from y that are missing in x
        for tra in missing_x:
            y_data = y_data.set_index('transformation')
            y_data = y_data.drop(tra, axis=0)
            y_data = y_data.reset_index()

    elif len(x) == len(y):
        x_data = x_data
        y_data = y_data
    
    eng_dict_plot_exp_for_eng[eng] = x_data
    eng_dict_plot_eng[eng] = y_data

    #make sure each x and y reflect the new data
    # x = (np.array(x_data.iloc[0:len(x_data),1])).reshape(-1,1)
    # y = np.array(y_data.iloc[0:len(y_data), 1])

    x = (np.array(x_data.iloc[0:len(x_data),1]))
    y = np.array(y_data.iloc[0:len(y_data), 1])

    #bootstrap to get measure of confidence intervals for these values
    #make a df of the values
    data_df = pd.DataFrame({'x': x, 'y': y})

    x = (data_df['x'].to_numpy()).reshape(-1,1)
    y = data_df['y'].to_numpy()

    #make the model and fit the data
    model = LinearRegression().fit(x,y)
    slope = model.coef_
    intercept = model.intercept_
    r_sq = model.score(x,y)
    mue = (abs(data_df['y'] - data_df['x']).sum())/len(x_data)

    # resample with replacement each row
    boot_slopes = []
    boot_interc = []
    boot_r_sq = [] 
    boot_mues = []
    #number of bootstrapping samples
    n_boots = 1000

    for n in range(n_boots):
        sample_df = data_df.sample(n=len(x_data), replace=True)
        # fit a linear regression
        x_temp = (sample_df['x'].to_numpy()).reshape(-1,1)
        y_temp = sample_df['y'].to_numpy()
        model_temp = LinearRegression().fit(x_temp,y_temp)
        slope_temp = model_temp.coef_
        intercept_temp = model_temp.intercept_
        r_sq_temp = model_temp.score(x_temp,y_temp) 

        mue_temp = (abs(sample_df['y'] - sample_df['x']).sum())/len(x_data)

        # append coefficients
        boot_interc.append(intercept_temp)
        boot_slopes.append(slope_temp)
        boot_r_sq.append(r_sq_temp)
        boot_mues.append(mue_temp)

    # #calculate the standard error for each
    # slope_err = (np.std(boot_slopes))/(math.sqrt(n_boots))
    # interc_err = (np.std(boot_interc))/(math.sqrt(n_boots))
    # mue_err = (np.std(boot_mues))/(math.sqrt(n_boots))
    # r_sq_err = (np.std(boot_r_sq))/(math.sqrt(n_boots))

    #calculate the standard deviation which tbh seems like thats what it was doing in the one from loeffler for these
    # deviation instead of the standard error
    slope_err = (np.std(boot_slopes))
    interc_err = (np.std(boot_interc))
    mue_err = (np.std(boot_mues))
    r_sq_err = (np.std(boot_r_sq))


    #append data frame 
    equations = equations.append({'engine' : eng ,'slope' : slope, 'slope_err' : slope_err, \
                                  'intercept' : intercept, 'intercept_err' : interc_err, \
                                    'r2' : r_sq, 'r2_err' : r_sq_err, \
                                    'MUE': mue, 'MUE_err': mue_err}, ignore_index=True).dropna(how='all')
    
    print(equations)
    #clear y data and missing_y incase so free for next engine. reset x data to the original.
    del y_data
    del missing_y
    x_data = eng_dict_val['experimental']

equations.to_csv(f'{res_dir}/linear_regression_equations.csv', index=False) #Write to output for later incase

#MAE between implementations
MAE_matrix = pd.DataFrame(columns=['AMBER','AMBER_err','SOMD','SOMD_err','GROMACS','GROMACS_err'], index=engine)

for eng in engine:
    x_data = eng_dict_val[eng]
    for eng2 in engine:
        y_data = eng_dict_val[eng2]

        x = np.array(x_data.iloc[0:len(x_data),1])
        y = np.array(y_data.iloc[0:len(y_data),1])

        data_df = pd.DataFrame({'x': x, 'y': y})

        x = (data_df['x'].to_numpy()).reshape(-1,1)
        y = data_df['y'].to_numpy()

        #mae = (abs(data_df['y'] - data_df['x']).sum())/len(x_data)
        mae = (abs(data_df['y'] - data_df['x']))
         
        #number of bootstrapping samples
        n_boots = 1000
        boot_maes = []
        for n in range(n_boots):
            sample_df = data_df.sample(n=len(x_data), replace=True)
            x_temp = (sample_df['x'].to_numpy()).reshape(-1,1)
            y_temp = sample_df['y'].to_numpy()
            mae_temp = []
            #mae_temp = (abs(sample_df['y'] - sample_df['x']).sum())/len(x_data)
            for i in range(len(mae)-1):
                r = random.randint(0,len(mae)-1)
                mae_temp.append(mae[r])
            mae_temp_mean = np.mean(mae_temp)
            # append coefficients
            boot_maes.append(mae_temp_mean)

        #calculate the standard error for each
        #mae_err = (np.std(boot_maes))/(math.sqrt(n_boots))
        #but this to calc std
        mae_err = (np.std(boot_maes))
        mae_mean = np.mean(boot_maes)

        MAE_matrix.loc[eng, eng2] = mae_mean
        MAE_matrix.loc[eng, f'{eng2}_err'] = mae_err

MAE_matrix.to_csv(f'{res_dir}/MAE_matrix.csv')


###

# plotting all together. Using the dictionaries saved before for plotting

plt.rc('font', size=12)
fig, ax = plt.subplots(figsize=(8,8))

colour = ['orange','orchid','darkturquoise','midnightblue']
lines = []

for eng,col in zip(engine,colour):
    x_data = eng_dict_plot_exp_for_eng[eng]
    y_data = eng_dict_plot_eng[eng]
    x = (np.array(x_data.iloc[0:len(x_data),1])).reshape(-1,1)
    y = np.array(y_data.iloc[0:len(y_data), 1])
    
    # TODO changed shapes here for the different proteins
    scatterplot = [plt.scatter(x[:4], y[:4], zorder=10, c=col, label="TYK2"),
                   plt.scatter(x[4:5], y[4:5], zorder=10, c=col, marker="D", label="p38"),
                   plt.scatter(x[5:], y[5:], zorder=10, c=col, marker="s",label="MCL1")]
    lines += plt.plot(0,0,c=col, label=eng)

    #plotting error bars
    y_er = np.array((eng_dict_plot_eng[eng]).iloc[0:len(eng_dict_plot_eng[eng]), 2])
    # x_er = np.array((eng_dict_plot_exp_for_eng[eng]).iloc[1:len(eng_dict_plot_exp_for_eng[eng]), 2])
    plt.errorbar(x , y,
                yerr=y_er,
                # xerr=x_er,   # comment this line to hide experimental error bars \
                            # as this can sometimes overcrowd the plot.
                ls="none",
                lw=0.5, 
                capsize=2,
                color="black",
                zorder=5
                )

    #plotting lines - need to change this part so incl from the correct written equations if want a linear fit line
    # x_line = np.linspace(-2,2,20)
    # y_line = (slope)*(x_line) + (intercept)
    # ax.plot(x_line, y_line, label=eng)

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, loc='upper left')
# plt.legend(scatterplot, ["TYK2","p38","MCL1"])

# plot 1/2 kcal bounds:
plt.fill_between(
                x=[-100, 100], 
                y2=[-100.25,99.75],
                y1=[-99.75, 100.25],
                lw=0, 
                zorder=-10,
                alpha=0.3,
                color="grey")
# upper bound:
plt.fill_between(
                x=[-100, 100], 
                y2=[-99.5,100.5],
                y1=[-99.75, 100.25],
                lw=0, 
                zorder=-10,
                color="grey", 
                alpha=0.2)
# lower bound:
plt.fill_between(
                x=[-100, 100], 
                y2=[-100.25,99.75],
                y1=[-100.5, 99.5],
                lw=0, 
                zorder=-10,
                color="grey", 
                alpha=0.2)

# get the bounds. This can be done with min/max or simply by hand.
all_freenrg_values_pre = []
for eng in engine:
    x_data = eng_dict_plot_exp_for_eng[eng]
    y_data = eng_dict_plot_eng[eng]
    x = (np.array(x_data.iloc[0:len(x_data),1])).tolist()
    y = (np.array(y_data.iloc[0:len(y_data), 1])).tolist()
    all_freenrg_values_pre.append(x)
    all_freenrg_values_pre.append(y)

all_freenrg_values = []
for sublist in all_freenrg_values_pre:
    for item in sublist:
        all_freenrg_values.append(item)

min_lim = min(all_freenrg_values)   
max_lim = max(all_freenrg_values)

# for a scatterplot we want the axis ranges to be the same. 
plt.xlim(min_lim*1.3, max_lim*1.3)
plt.ylim(min_lim*1.3, max_lim*1.3)

plt.axhline(color="black", zorder=1)
plt.axvline(color="black", zorder=1)

#plt.xlabel('ΔΔG for experimental (kcal/mol)')
#plt.ylabel('ΔΔG for calculated (kcal/mol)')
plt.ylabel("Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
plt.xlabel("Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")

#plt.title('calculated and experimental linear regression')

plt.savefig(f'{res_dir}/r2_correlation.png')






###

#can currently do in excel by opening the csv file of results

#plotting a bar graph

#tbh best to make new data frames so that they're also together with the correct transformation so in the data frame trans exp each eng
#also want exp from but need to remove all the empty ones if there are any as this has not been done yet - do I need to? actually no 

# exper_bar = eng_dict_val['experimental']

# try:
#     exper_bar = exper_bar.set_index('freenrg')
#     try:
#         exper_bar = exper_bar.drop("NaN", axis=0)
#     except:
#         pass
#     exper_bar = exper_bar.reset_index()
#     exper_bar = exper_bar.set_index('transformation')
#     exper_bar = exper_bar.reset_index()
# except:
#     pass

# bar_data = pd.DataFrame(index=["transformation", "experimental", "exp_err", "SOMD", "SOMD_err", "AMBER", "AMBER_err", "GROMACS", "GROMACS_err"])

# for tra in trans:

#     #df = exper_bar
#     df = eng_dict_val['experimental']
#     try:
#         if df.loc[df['transformation'] == tra, 'transformation'].iloc[0] == tra :
#             exp_freenrg = df.loc[df['transformation'] == tra, 'freenrg'].iloc[0]
#             exp_error = df.loc[df['transformation'] == tra, 'error'].iloc[0]
#         else:
#             exp_freenrg = 0
#             exp_error = 0
#     except:
#         exp_freenrg = 0
#         exp_error = 0

#     try:
#         df = eng_dict_val['SOMD']
#         if df.loc[df['transformation'] == tra, 'transformation'].iloc[0] == tra :
#             somd_freenrg = df.loc[df['transformation'] == tra, 'freenrg'].iloc[0]
#             somd_error = df.loc[df['transformation'] == tra, 'error'].iloc[0]
#         else:
#             somd_freenrg = 0
#             somd_error = 0
#     except:
#         somd_freenrg = 0
#         somd_error = 0

#     try:
#         df = eng_dict_val['AMBER']
#         if df.loc[df['transformation'] == tra, 'transformation'].iloc[0] == tra :
#             amber_freenrg = df.loc[df['transformation'] == tra, 'freenrg'].iloc[0]
#             amber_error = df.loc[df['transformation'] == tra, 'error'].iloc[0]
#         else:
#             amber_freenrg = 0
#             amber_error = 0
#     except:
#         amber_freenrg = 0
#         amber_error = 0

#     try:
#         df = eng_dict_val['GROMACS']
#         if df.loc[df['transformation'] == tra, 'transformation'].iloc[0] == tra :
#             gromacs_freenrg = df.loc[df['transformation'] == tra, 'freenrg'].iloc[0]
#             gromacs_error = df.loc[df['transformation'] == tra, 'error'].iloc[0]
#         else:
#             gromacs_freenrg = 0
#             gromacs_error = 0
#     except:
#         gromacs_freenrg = 0
#         gromacs_error = 0
    
#     bar_data = bar_data.append({"transformation":tra,\
#                                  "SOMD":somd_freenrg, "SOMD_err":somd_error,\
#                                  "AMBER":amber_freenrg, "AMBER_err":amber_error,\
#                                  "GROMACS":gromacs_freenrg, "GROMACS_err":gromacs_error,\
#                                  "experimental":exp_freenrg, "exp_err":exp_error}, ignore_index=True).dropna(how='all')

# bar_data.to_csv('all_results_table.csv', index=False) #can write to csv so all together tooo

# #np.nan need to incooperate somehow earlier
# #bar is broken in the sense that it currently only works 

# #need to substitute the NaN as well for 0 - because thats how they will be written in the table and cant plot those
# #bar_data.replace(to_replace='NaN', value=0) there must be a better way to do it but will simply open the file and replace
# fin = open('all_results_table.csv', 'r')
# fiout = open('all_results_table_bar.csv', 'w')
# for line in fin:
#     fiout.write(line.replace(',,',',0.0,'))
# fin.close()
# fiout.close()
# ##this doesnt currently work at all, manually added 0 to the all results table bar
# #bar_data = pd.read_csv('all_results_table_bar.csv',sep=',', header='infer', index_col=None)
# # print(bar_data['SOMD'])
# # print(bar_data.loc['SOMD'].index=False)

# bar_data = pd.read_csv('all_results_table_bar.csv',sep=',', header='infer', index_col='transformation')

# import matplotlib.pyplot as plot

# bar_data.plot.bar()

# plot.savefig('bar_all.png')

# #actually plotting
# width = 0.15
# offset = [(- width*1.5),(- width*0.5),(+ width*0,5),(+ width*1.5)]

# fig, ax = plt.subplots(figsize=(10,10))

# #better way to append exp to engines used but thats currently not workign so
# engine_bar = ['SOMD', 'AMBER', 'GROMACS', 'experimental']
# #engine_bar = ['AMBER', 'experimental']
# # engine_bar = engine_bar.append('experimental')
# # print(f'engine ins {engine}')

# error_dict = {'AMBER' : bar_data['AMBER_err'], 'SOMD' : bar_data['SOMD_err'],\
#               'GROMACS' : bar_data['GROMACS_err'], 'experimental' : bar_data['exp_err']}
# print(error_dict)
# x_locs = np.arange(len(trans))

# # for eng,col,off in zip(engine_bar,colour,offset):
# #     y_er = error_dict[eng]
# #     if eng == 'experimental':
# #         ax.bar(x_locs + off, bar_data[eng], width=width, label=eng, color=col)
# #     else:
# #         ax.bar(x_locs + off, bar_data[eng], width=width, yerr= y_er, label=eng, color=col)
        
# offset = [(- width*1.5),(- width*0.5),(+ width*0,5),(+ width*1.5)]
# colour = ['orange','orchid','darkturquoise','midnightblue']

# ax.bar(x_locs + (width*1.5), bar_data['experimental'], width=width, label='experimental', color='midnightblue')
# #ax.bar(x_locs + offset[1], bar_data['AMBER'], width=width, yerr= error_dict['AMBER'], label='AMBER', color=col[1])
        

# plt.axhline(color="black")
# plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
# plt.xticks(x_locs, trans, rotation=50, ha="left")
# plt.legend()

# plt.savefig('bar_all.png')

# #should be according to the maximum length ie no of transformations
# x_locs = np.arange(len(trans))

# for eng,col,off in zip(engine,colour,offset):
#     y = (eng_dict[eng]).iloc[0:len(eng_dict_val[eng]), 2]
#     if eng == 'experimental':
#         ax.bar(x_locs + off, y, width=width, label=eng, color=col)
#     else:
#         ax.bar(x_locs + off, y, width=width, yerr=y_er, label=eng, color=col)

# plt.axhline(color="black")
# plt.ylabel("$\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
# plt.xticks(x_locs, trans, rotation=50, ha="right")
# plt.legend()

# plt.savefig('bar_all.jpg')

# ###


os.chdir(res_dir)