#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

"""
@author: Dilek Koptekin
"""

import pandas as pd
import numpy as np
from itertools import product
import os
import argparse

def create_combinations(model, n_pops):
    if n_pops == 3:
        return pd.DataFrame(product(targets,sources[model[1]], sources[model[2]]))
    if n_pops == 4:
        return pd.DataFrame(product(targets,sources[model[1]], sources[model[2]], sources[model[3]]))
    if n_pops >= 5:
        return pd.DataFrame(product(targets,sources[model[1]], sources[model[2]], sources[model[3]], sources[model[4]]))

def create_parfile(left_name, par_name, par_lines, outgroup):
    with open(par_name, "w") as par:
        for line in parfile:
            if "popleft" in line or "popright" in line:
                ll = line.split()
                if ll[0].startswith("popleft"):
                    ll[1] = left_name
                if ll[0].startswith("popright"):
                    ll[1] = outgroup
                line = " ".join(ll) + "\n"
            par.write(line)

def create_sbatch(run_list, host_name, n_threads, outputdir, admixtools_path, no):
    sbatch_param = f"#!/bin/bash -l\n#SBATCH -J qpAdm \n#SBATCH -p {host_name}\n#SBATCH -n {n_threads}\n#SBATCH -t 5-00:00:00\n#SBATCH -D {outputdir} \n#SBATCH -o slurm-%j-%N-%u.out\n#SBATCH -e slurm-%J-%N-%u.err\nexport PATH=$PATH:{admixtools_path}\n"
    qpAdm_run = "\n".join(run_list)
    sbatch_name = "sbatch_qpAdm" + "_" + info + "_" + host_name + "_t" + str(n_threads) +"_"+ str(no).zfill(4) + ".sh"
    sbatch_out = sbatch_param + qpAdm_run + "\nwait"
    with open(sbatch_name, "w") as sbatch:
        sbatch.write(sbatch_out)

def create_bash(run_list, n_threads, outputdir, admixtools_path, no):
    bash_param = f"#!/bin/bash\ncd  {outputdir}  \nexport PATH=$PATH:{admixtools_path}\n"
    qpAdm_run = "\n".join(run_list)
    bash_name = "qpAdm" + "_" + info + "_t" + str(n_threads) +"_"+ str(no).zfill(4) + ".sh"
    bash_out = bash_param + qpAdm_run
    with open(bash_name, "w") as bash:
        bash.write(bash_out)

def create_result_df(n_col, n_row):
    colnames = ['target', 'pop1', 'pop1_mixture', 's.e_1', 'z_pop1','pop2', 'pop2_mixture', 's.e_2', 'z_pop2','pop3', 'pop3_mixture', 's.e_3', 'z_pop3', 'pop4', 'pop4_mixture', 's.e_4', 'z_pop4', 'p_value', 'outgroup', 'feasible', 'z_eval',  'logfile_path']
    results = pd.DataFrame(columns = colnames, index = range(n_row))
    str_col = ['target', 'pop1', 'pop2', 'pop3', 'pop4', 'outgroup','feasible','z_eval', 'logfile_path']
    int_col = results.columns.difference(str_col)
    results[int_col] = results[int_col].apply(pd.to_numeric)
    return results

def eval_z(row):
    z = ['z_pop1', 'z_pop2', 'z_pop3', 'z_pop4']
    z = row[z].dropna()
    zval = []
    for i in range(len(z)):
        if z[i] < 2:
            zval.append("fail")
        else:
            zval.append("pass")
    if "fail" in zval:
        return "fail"
    else:
        return "pass"

def eval_p(row):
    if (row['p_value'] < 0.01):
        return 'p-value/not feasible'
    elif (row['p_value'] > 0.01) & (row['z_eval'] == "fail"):
        return 'propotions/Z<2'
    elif (row['p_value'] > 0.01) & (row['z_eval'] == "pass"):
        return 'feasible'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser=argparse.ArgumentParser(description='Running qpAdm')
    parser.add_argument("-i","--input", type=str,  dest='info', required=True, help="")
    parser.add_argument("-o","--outgroup_name", type=str,  dest='outgroup_name', help="")
    parser.add_argument("-p","--par", type=str,  dest='default_par', required=True, help="")
    parser.add_argument("-w","--write", default=True, type=bool)
    parser.add_argument("-r","--read", default=False, type=bool)
    parser.add_argument("--run_mode", dest='mode', type=str,  choices=['sbatch', 'bash'], default="sbatch",help="")
    parser.add_argument("--admixtools_path", dest='admixtools_path', default="/usr/local/sw/AdmixTools-7.0.2/bin/", type=str, help="")
    parser.add_argument("-t","--threads", dest='threads', default=1, type=int, help="")
    parser.add_argument("-ho","--host_name", dest='host', default="chimp", type=str, help="")
    parser.add_argument("--out_path", dest='output_path',  default="./", type=str)

    args = parser.parse_args()
    info = args.info
    default_par = args.default_par

    # for sbatch file
    n_threads = args.threads
    host=args.host
    outdir=args.output_path
    adm_path = args.admixtools_path
    #adm_path = "/usr/local/sw/AdmixTools-7.0.2/bin/"
    os.chdir(outdir)
    # read pop-model file
    df = pd.read_table(info,
                       header = 0, index_col = 0,
                       sep = '\t').transpose()

    names = df.index
    df_models = pd.DataFrame(df["model"])
    df_pops = pd.DataFrame(df["pop"])
    n_models = df_models.shape[1]
    models = {i+1:list(df_models.iloc[:,i].dropna()) for i in range(n_models)}
    targets = list(df["pop"].iloc[0].dropna())
    sources = {names[i]:list(df["pop"].iloc[i].dropna()) for i in range(1,len(names))}


    parfile = open(default_par, "r").readlines()

    if args.outgroup_name is not None:
        outgroup=args.outgroup_name
    else:
        outgroup=parfile[4].split()[1]

    ind_filename = parfile[2].split()[1]
    ind_pops = np.unique(np.loadtxt(ind_filename, usecols=(2,), dtype = str))
    right_pops = np.loadtxt(outgroup, usecols=(0,), dtype = str)
    model_pops = np.append(targets, right_pops)
    for source_list in sources.values():
        model_pops = np.append(source_list, model_pops)
    no_ind = []
    for model_pop in model_pops:
        if model_pop not in ind_pops:
            no_ind.append(model_pop)
    if len(no_ind) > 0:
        print(f"The populations below were not found in the ind file")
        print("\n".join(no_ind))

    if args.read:
        args.write = False

    if args.write:
        # prepare left and par files
        par_list = []
        qpAdm_list_all = []
        poplist = []
        for i in models:
            model = models[i]
            n_pops = len(model)
            comb = create_combinations(model, n_pops)
            for run in range(comb.shape[0]):
                filename_left = (info + '_model' + str(i) + "_left_" + str(run+1).zfill(4))
                left = comb.iloc[run]
                if pd.unique(left).shape == left.shape:
                    poplist.append("\t".join(list(left)))
                    left.to_csv(filename_left, index=False, header = False, sep = "\n")
                    filename_par = "par_qpAdm_"+ outgroup + "_" + "_".join(list(left))
                    filename_log = "log_qpAdm_" + outgroup + "_" + "_".join(list(left))
                    create_parfile(filename_left, filename_par, parfile, outgroup)
                    par_list.append(filename_par)
                    qpAdm_list_all.append("(qpAdm -p " + filename_par + " > " + filename_log + ") &")
                else:
                    print("The following model includes duplicate population in left \n ")
                    print(" ".join(left))
        if args.mode == "sbatch":
            if len(qpAdm_list_all) > n_threads:
                [create_sbatch(qpAdm_list_all[no_run:no_run+n_threads],host, n_threads, outdir, adm_path, no_run) for no_run in range(0, len(qpAdm_list_all), n_threads)]
            else:
                create_sbatch(qpAdm_list_all, host, n_threads, outdir, adm_path, 0)
        if args.mode == "bash":
            if len(qpAdm_list_all) > n_threads:
                [create_bash(qpAdm_list_all[no_run:no_run+n_threads], n_threads, outdir, adm_path, no_run) for no_run in range(0, len(qpAdm_list_all), n_threads)]
            else:
                create_bash(qpAdm_list_all, n_threads, outdir, adm_path, 0)
        parlist_filename = "all_run_" + info + "_" + host + "_t" + str(n_threads)
        poplist_filename = "all_combination_" + info
        with open(parlist_filename, "w") as par_out:
            par_out.write("\n".join(par_list))
        with open(poplist_filename, "w") as poplist_out:
            poplist_out.write("\n".join(poplist))

    elif args.read:
        run_comb = open(("all_combination_" + info), "r").readlines()
        results = create_result_df(len(names), len(run_comb))
        results['outgroup'] = outgroup
        for l,log in enumerate(run_comb):
            log_pops = log.split()
            log_name = "log_qpAdm_" + outgroup + "_" + "_".join(list(log_pops))
            if os.path.isfile(log_name):
                results.at[l,'logfile_path'] = os.path.join(os.getcwd(), log_name)
                log_file = open(log_name,'r').readlines()
                for no, res in enumerate(log_file):
                    if 'left pops:' in res:
                        results.at[l,'target'] = log_file[no+1].rstrip('\n')
                        results.at[l,'pop1'] = log_file[no+2].rstrip('\n')
                        results.at[l,'pop2'] = log_file[no+3].rstrip('\n')
                        if len(log_pops) >= 4:
                            results.at[l,'pop3'] = log_file[no+4].rstrip('\n')
                            if len(log_pops) >= 5:
                                results.at[l,'pop4'] = log_file[no+5].rstrip('\n')
                    if 'best coefficients:' in res:
                        percent = res.split()[2:]
                        results.at[l,'pop1_mixture'] = float(percent[0])
                        results.at[l,'pop2_mixture'] = float(percent[1])
                        if len(log_pops) >= 4:
                            results.at[l,'pop3_mixture'] = float(percent[2])
                            if len(log_pops) >= 5:
                                results.at[l,'pop4_mixture'] = float(percent[3])
                    if 'std. errors:' in res:
                        std = res.split()[2:]
                        results.at[l,'s.e_1'] = float(std[0])
                        results.at[l,'s.e_2'] = float(std[1])
                        if len(log_pops) >= 4:
                            results.at[l,'s.e_3'] = float(std[2])
                            if len(log_pops) >= 5:
                                results.at[l,'s.e_4'] = float(std[3])
                    if 'full rank' in res:
                        results.at[l,'p_value'] = float(log_file[no+2].split()[-1])


        results['z_pop1'] = results['pop1_mixture'] / results['s.e_1']
        results['z_pop2'] = results['pop2_mixture'] / results['s.e_2']
        results['z_pop3'] = results['pop3_mixture'] / results['s.e_3']
        results['z_pop4'] = results['pop4_mixture'] / results['s.e_4']
        results['z_eval'] = results.apply(eval_z, axis=1).astype(str)
        results['feasible'] = results.apply(eval_p, axis=1).astype(str)
        results.to_csv(("results_poplist_" + info + ".txt"), sep ='\t', index = False)
