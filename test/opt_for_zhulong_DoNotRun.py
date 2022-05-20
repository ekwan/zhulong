# Modification of last section "run the experiments" in zhulong.py:

parameter_bounds = parameter_space.get_bounds_dict()
method = 'BO' # or 'RS' or 'LM'
initMethod = 'batch8' # 'random' or 'batch8' for LM and BO
n_total = 24 # total iterations

dat = pd.DataFrame()
f_best = 0
for r in range(n_total):
    # get a set of new parameters for round r
    para_r = optimization(dat, parameter_bounds,method, seed = 1234, initMethod = initMethod) 
    df_para_r = pd.DataFrame({k: [v] for k, v in para_r.items()})
    # define an experiment object for this round
    experiment_r = Experiment.create(
                         parameter_space,                   
                         solvent=para_r['Solvent'],                   
                         temperature=para_r['Temperature'],                    
                         starting_material_volume=100,     
                         reagent=para_r['Reagent'],                   
                         reagent_equivalents=para_r['Reagent_equiv'],        
                         additive=para_r['Additive'], 
                         additive_mole_percent=para_r['AdditiveLoading'], 
                         light_stage=para_r['Stage']) 
    print("running experiment")
    print(experiment_r)
    chemspeed.run_experiment(experiment_r)
    history = experiment_r.history
    print("final result is:")
    print(f"{history['times']=}")
    print(f"{history['values']=}")
    f_r = history['values'][-1]
    f_best = max(f_best, f_r)
    df_para_r = df_para_r.assign(f = f_r, f_best = f_best)
    dat = pd.concat([dat,df_para_r],ignore_index = True)
    print("Round: ", r+1, "f: ", f_r,"f_best",f_best,"\n")

dat = dat.assign(Round = range(1,n_total+1))
dat = dat.assign(method = method)
dat.to_csv("dat.csv")
