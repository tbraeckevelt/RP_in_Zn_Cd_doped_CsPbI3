from pathlib import Path
import parsl
from parsl.data_provider.files import File
import psiflow
import os
import molmod.units
import numpy as np
from ownscripts.bash_app_python import bash_app_python
from ownscripts.utils import run_MD

from utils_debug import plotvol, writeatoms, printenergy, writeinputtraj, create_ave_struc, get_bondangles_distribution

app_plotvol = bash_app_python(plotvol, precommand = "", executors=['default_threads'])
app_writeatoms = bash_app_python(writeatoms, precommand = "", executors=['default_threads'])
app_printenergy = bash_app_python(printenergy, precommand = "", executors=['default_threads'])
app_writeinputtraj = bash_app_python(writeinputtraj, precommand = "", executors=['default_threads'])
app_create_ave_struc = bash_app_python(create_ave_struc, precommand = "", executors=['default_threads'])
app_get_bondangles_distribution = bash_app_python(get_bondangles_distribution, precommand = "", executors=['default_threads'])
app_run_MD = bash_app_python(run_MD, precommand = "", executors=['ModelEvaluation'])

def app_run_MD_general(output_folder, simname, path_atoms, mult_steps = 1, seed = 0):
    #path calculator
    path_calc = File(str(data_folder / "macemodels" / "OwnMLP_MACE.model")) # The model is currently compressed via gzip due to memory constraints of GitHub,
                                                                            # you can alos use then universal foundation model, mace mp 0
    #General settings MD
    total_steps = 100000 # for some simulation this number is doubled
    freq_step = 100
    temperature = 300
    pressure = 0.1*molmod.units.pascal*(10**6)
    
    path_output = File(str(output_folder / str("traj_"+simname+".xyz")))
    if os.path.exists(path_output.filepath) == False:
        MDapp = app_run_MD(
            execution_folder = output_folder,
            stderr           = str(output_folder / str("error_MD_"+simname+".txt")),
            stdout           = str(output_folder / str("output_MD_"+simname+".txt")),
            inputs           = [path_atoms, path_calc], 
            outputs          = [path_output], 
            steps            = int(mult_steps*total_steps), 
            step             = freq_step, 
            temperature      = temperature, 
            pressure         = pressure, 
            seed             = seed, 
            device           = "cuda")
        path_output = MDapp.outputs[0]
    return path_output

def app_create_ave_struc_general(output_folder, path_traj, simname):
    #General settings MD
    calib = 100 # number of printed steps (so multiply with freq_step for total steps) used as calibration and therefor not included in any post analysis
    path_ave_atoms = File(str(output_folder / str("ave_atoms_"+simname+".cif")))
    if os.path.exists(path_ave_atoms.filepath) == False:
        ave_at_app = app_create_ave_struc(
            execution_folder = output_folder,
            stderr           = str(output_folder / str("error_aveat_"+simname+".txt")),
            stdout           = str(output_folder / str("output_aveat_"+simname+".txt")),
            inputs           = [path_traj],
            outputs          = [path_ave_atoms],
            calib            = calib)
        return ave_at_app

def app_writeatoms_general(output_folder, path_traj, simname, index = -1):
    if index == -1:
        simname = "end_" + simname
    elif index == 0:
        simname = "begin_" + simname
    else:
        simname = "at"+str(index)+"_" + simname
    path_output = File(str(output_folder / str("atoms_"+simname+".xyz")))
    if os.path.exists(path_output.filepath) == False:
        wa_app = app_writeatoms(
            execution_folder = output_folder,
            stderr           = str(output_folder / str("error_write_"+simname+".txt")),
            stdout           = str(output_folder / str("output_write_"+simname+".txt")),
            inputs           = [path_traj], 
            outputs          = [path_output],
            index            = index)
        path_output = wa_app.outputs[0]
    return path_output


if __name__ == '__main__':
    psiflow.load()

    create_input_struc = True
    restart = True
    num_res = 1
    new_input_traj=[]
    ave_at_apps_lst = []
    calib = 100 # number of printed steps (so multiply with freq_step for total steps) used as calibration and therefor not included in any post analysis

    #Create input and output folders
    Path_folder = Path.cwd() 
    data_folder = Path_folder / "data"
    main_output_folder = Path_folder / "output"
    main_output_folder.mkdir(exist_ok=True)
    print("Start constructing the simulations")

    #Perovsktite simulations
    pv_data_folder = data_folder / "PerovskitePhase"
    pv_output_folder = main_output_folder / "PerovskitePhase"
    pv_output_folder.mkdir(exist_ok=True)
    print("submit perovskite simualtions")
    
    el_lst = ["Cd", "Zn"]
    subst_lst = ["Cs", "Pb"]
    ns_lst = range(9)
    sd_lst = range(10)

    #MD on pure perovskite phase - start alpha
    path_atoms = File(str(pv_data_folder / str("atoms_0.xyz")))
    pure_pv_outputtraj_dct = {}
    for sd in sd_lst:
        simname = "0_"+str(sd)
        pure_pv_outputtraj_dct[sd] = app_run_MD_general(pv_output_folder, simname, path_atoms, seed = sd)
        ave_at_apps_lst.append(app_create_ave_struc_general(pv_output_folder, pure_pv_outputtraj_dct[sd], simname))

    #MD on pure perovskite phase - start gamma
    path_atoms = File(str(pv_data_folder / str("atoms_gamma_0.xyz")))
    simname = "gamma_0"
    puregamma_pv_outputtraj = app_run_MD_general(pv_output_folder, simname, path_atoms)
    ave_at_apps_lst.append(app_create_ave_struc_general(pv_output_folder, puregamma_pv_outputtraj, simname))

    pv_outputtraj_dct = {}
    #Doped perovskites
    for el in el_lst:
        pv_outputtraj_dct[el]={}
        for subst in subst_lst:
            pv_outputtraj_dct[el][subst]={}
            for ns in ns_lst:
                pv_outputtraj_dct[el][subst][ns]={}
                path_atoms = File(str(pv_data_folder / str("atoms_"+el+"_"+subst+"_"+str(ns)+".xyz")))
                for sd in sd_lst:
                    if ns == 0:
                        pv_outputtraj_dct[el][subst][ns][sd] = pure_pv_outputtraj_dct[sd]
                    else:
                        simname = el+"_"+subst+"_"+str(ns)+"_"+str(sd)
                        pv_outputtraj_dct[el][subst][ns][sd] = app_run_MD_general(pv_output_folder, simname, path_atoms, seed = sd)
                        ave_at_apps_lst.append(app_create_ave_struc_general(pv_output_folder, pv_outputtraj_dct[el][subst][ns][sd], simname))
                    if create_input_struc and ns >= 1 and sd == 0:
                        new_input_traj.append(app_writeatoms_general(pv_output_folder, pv_outputtraj_dct[el][subst][ns][sd], simname))

    #Plot volume expansion:
    inputs = []
    for el in el_lst:
        for subst in subst_lst:
            for ns in ns_lst:
                for sd in sd_lst:
                    inputs.append(pv_outputtraj_dct[el][subst][ns][sd])
    volplotfile = File(str(pv_output_folder / str("Vol_plot.pdf")))
    if os.path.exists(volplotfile.filepath) == False:
        pv_plotvol_app = app_plotvol(
            execution_folder = pv_output_folder,
            stderr           = str(pv_output_folder / str("error_plotvol.txt")),
            stdout           = str(pv_output_folder / str("output_plotvol.txt")),
            inputs           = inputs, 
            outputs          = [volplotfile], 
            el_lst           = el_lst, 
            subst_lst        = subst_lst, 
            ns_lst           = ns_lst, 
            sd_lst           = sd_lst, 
            calib            = calib)

    if restart:
        #MD on pure perovskite phase - start alpha
        pure_pv_outputtraj_dct_res = {}
        for sd in sd_lst:
            simname = "0_"+str(sd) + "_res"
            pure_pv_outputtraj_dct_res[sd] = app_run_MD_general(pv_output_folder, simname, pure_pv_outputtraj_dct[sd], seed = sd)
            ave_at_apps_lst.append(app_create_ave_struc_general(pv_output_folder, pure_pv_outputtraj_dct_res[sd], simname))

        pv_outputtraj_dct_res = {}
        #Doped perovskites
        for el in el_lst:
            pv_outputtraj_dct_res[el]={}
            for subst in subst_lst:
                pv_outputtraj_dct_res[el][subst]={}
                for ns in ns_lst:
                    pv_outputtraj_dct_res[el][subst][ns]={}
                    for sd in sd_lst:
                        if ns == 0:
                            pv_outputtraj_dct_res[el][subst][ns][sd] = pure_pv_outputtraj_dct_res[sd]
                        else:
                            simname = el+"_"+subst+"_"+str(ns)+"_"+str(sd) + "_res"
                            pv_outputtraj_dct_res[el][subst][ns][sd]=app_run_MD_general(pv_output_folder,simname,pv_outputtraj_dct[el][subst][ns][sd],seed=sd)
                            ave_at_apps_lst.append(app_create_ave_struc_general(pv_output_folder, pv_outputtraj_dct_res[el][subst][ns][sd], simname))

        #Plot volume expansion:
        inputs = []
        for el in el_lst:
            for subst in subst_lst:
                for ns in ns_lst:
                    for sd in sd_lst:
                        inputs.append(pv_outputtraj_dct_res[el][subst][ns][sd])
        volplotfile = File(str(pv_output_folder / str("Vol_plot_res.pdf")))
        if os.path.exists(volplotfile.filepath) == False:
            pv_plotvol_app_res = app_plotvol(
                execution_folder = pv_output_folder,
                stderr           = str(pv_output_folder / str("error_plotvol_res.txt")),
                stdout           = str(pv_output_folder / str("output_plotvol_res.txt")),
                inputs           = inputs, 
                outputs          = [volplotfile], 
                el_lst           = el_lst, 
                subst_lst        = subst_lst, 
                ns_lst           = ns_lst, 
                sd_lst           = sd_lst, 
                calib            = calib)


    #RP simulations
    rp_data_folder = data_folder / "RPphase"
    rp_output_folder = main_output_folder / "RPphase"
    rp_output_folder.mkdir(exist_ok=True)
    print("submit RP simualtions")
    rp_outputtraj_dct ={}

    #Cs substitution
    rp_outputtraj_dct["Cs"] = {}
    plane_lst = ["001", "110"]
    el_lst = ["Cd", "Zn"]
    subst_lst = range(6)
    sd_lst = range(4)
    for plane in plane_lst:
        rp_outputtraj_dct["Cs"][plane] = {}
        for el in el_lst:
            rp_outputtraj_dct["Cs"][plane][el] = {}
            for subst in subst_lst:
                rp_outputtraj_dct["Cs"][plane][el][subst] = {}
                path_atoms = File(str(rp_data_folder / str("atoms_"+plane+"_"+el+"_Cs"+str(subst)+".xyz")))
                for sd in sd_lst:
                    simname = plane+"_"+el+"_Cs"+str(subst)+"_"+str(sd)
                    rp_outputtraj_dct["Cs"][plane][el][subst][sd] = app_run_MD_general(rp_output_folder, simname, path_atoms, seed = sd)
                    ave_at_apps_lst.append(app_create_ave_struc_general(rp_output_folder, rp_outputtraj_dct["Cs"][plane][el][subst][sd], simname))
                    if create_input_struc and sd == 0:
                        new_input_traj.append(app_writeatoms_general(rp_output_folder, rp_outputtraj_dct["Cs"][plane][el][subst][sd], simname))
                
    #Get average energy
    inputs = []
    for plane in plane_lst:
        for el in el_lst:
            for subst in subst_lst:
                for sd in sd_lst:
                    inputs.append(rp_outputtraj_dct["Cs"][plane][el][subst][sd])
    energy_txtfile = File(str(rp_output_folder / str("Summary_energy_Cspos.txt")))
    if os.path.exists(energy_txtfile.filepath) == False:
        rp_Cs_printenergy_app = app_printenergy(
            execution_folder = rp_output_folder,
            stderr           = str(rp_output_folder / str("error_printenergy_Cs.txt")),
            stdout           = str(rp_output_folder / str("output_printenergy_Cs.txt")),
            inputs           = inputs, 
            outputs          = [energy_txtfile], 
            plane_lst        = plane_lst, 
            el_lst           = el_lst, 
            subst_lst        = subst_lst, 
            sd_lst           = sd_lst, 
            calib            = calib)


    #Pb substitution
    rp_outputtraj_dct["Pb"] = {}
    el_lst_Pb = el_lst + ["Pb"]
    coord_lst = ["Pb4co", "Pb5co", "Pb6co"]
    sd_lst = range(4)
    for plane in plane_lst:
        rp_outputtraj_dct["Pb"][plane] = {}
        for el in el_lst_Pb:
            rp_outputtraj_dct["Pb"][plane][el] = {}
            for co in coord_lst:
                rp_outputtraj_dct["Pb"][plane][el][co] = {}
                path_atoms = File(str(rp_data_folder / str("atoms_"+plane+"_"+el+"_"+co+".xyz")))
                for sd in range(10):
                    simname = plane+"_"+el+"_"+co+"_"+str(sd)
                    rp_outputtraj_dct["Pb"][plane][el][co][sd] = app_run_MD_general(rp_output_folder, simname, path_atoms, seed = sd)
                    ave_at_apps_lst.append(app_create_ave_struc_general(rp_output_folder, rp_outputtraj_dct["Pb"][plane][el][co][sd], simname))
                    if create_input_struc and sd == 0:
                        new_input_traj.append(app_writeatoms_general(rp_output_folder, rp_outputtraj_dct["Pb"][plane][el][co][sd], simname))

    #Get average energy
    inputs = []
    for plane in plane_lst:
        for el in el_lst_Pb:
            for coord in coord_lst:
                for sd in sd_lst:
                    inputs.append(rp_outputtraj_dct["Pb"][plane][el][coord][sd])
    energy_txtfile = File(str(rp_output_folder / str("Summary_energy_Pbpos.txt")))
    if os.path.exists(energy_txtfile.filepath) == False:
        rp_Pb_printenergy_app = app_printenergy(
            execution_folder = rp_output_folder,
            stderr           = str(rp_output_folder / str("error_printenergy_Pb.txt")),
            stdout           = str(rp_output_folder / str("output_printenergy_Pb.txt")),
            inputs           = inputs, 
            outputs          = [energy_txtfile], 
            plane_lst        = plane_lst, 
            el_lst           = el_lst_Pb, 
            subst_lst        = coord_lst, 
            sd_lst           = sd_lst, 
            calib            = calib)

    
    #interface simulations
    if_data_folder = data_folder / "interface"
    if_output_folder = main_output_folder / "interface"
    if_output_folder.mkdir(exist_ok=True)
    print("submit interface simualtions")
    #7 different input struc - do MD:
    sim_lst = ["001RP_small", "001RP_uniform", "001RP_varying", "110RP_verysmall", "110RP_small", "110RP_normal", "110RP_big"]
    forinput_lst = ["001RP_small", "001RP_varying", "110RP_verysmall", "110RP_small"]
    if_outputtraj_dct= {}
    for simname in sim_lst:
        path_atoms = File(str(if_data_folder / str("atoms_"+simname+".xyz")))
        if_outputtraj_dct[simname] = app_run_MD_general(if_output_folder, simname, path_atoms, mult_steps = 2)
        ave_at_apps_lst.append(app_create_ave_struc_general(if_output_folder, if_outputtraj_dct[simname], simname))
        if create_input_struc and simname in forinput_lst:
            new_input_traj.append(app_writeatoms_general(if_output_folder, if_outputtraj_dct[simname], simname))

    #big and varying structures, that will need a restart
    simname = "110RP_big_varying"
    path_atoms = File(str(if_data_folder / str("atoms_110RP_big_varying.xyz")))
    if_outputtraj_dct[simname] = app_run_MD_general(if_output_folder, simname, path_atoms)
    ave_at_apps_lst.append(app_create_ave_struc_general(if_output_folder, if_outputtraj_dct[simname], simname))
    if restart:
        for nr in range(num_res):
            path_atoms = if_outputtraj_dct[simname]
            simname = "110RP_big_varying_res" + str(nr)
            if_outputtraj_dct[simname] = app_run_MD_general(if_output_folder, simname, path_atoms)
            ave_at_apps_lst.append(app_create_ave_struc_general(if_output_folder, if_outputtraj_dct[simname], simname))


    #Get Pb-I-Pb bond angles distribution
    bondangle_plot_file = File(str(main_output_folder / str("Bondangles_comparison.pdf")))
    labels = ["pure_perovskite", "pure_RP001", "pure_RP110", "interface_001_uniform", "interface_001_varying", "interface_110_big", "110RP_big_varying"]
    if restart:
        inputs = [pure_pv_outputtraj_dct_res[0], 
                rp_outputtraj_dct["Pb"]["001"]["Pb"]["Pb6co"][0],
                rp_outputtraj_dct["Pb"]["110"]["Pb"]["Pb6co"][0],
                if_outputtraj_dct["001RP_uniform"],
                if_outputtraj_dct["001RP_varying"],
                if_outputtraj_dct["110RP_big"],
                if_outputtraj_dct["110RP_big_varying_res" + str(num_res - 1)]]
    else:
        inputs = [pure_pv_outputtraj_dct[0], 
                rp_outputtraj_dct["Pb"]["001"]["Pb"]["Pb6co"][0],
                rp_outputtraj_dct["Pb"]["110"]["Pb"]["Pb6co"][0],
                if_outputtraj_dct["001RP_uniform"],
                if_outputtraj_dct["001RP_varying"],
                if_outputtraj_dct["110RP_big"],
                if_outputtraj_dct["110RP_big_varying"]]
    if os.path.exists(bondangle_plot_file.filepath) == False:
        ba_app = app_get_bondangles_distribution(
            execution_folder = main_output_folder,
            stderr           = str(main_output_folder / str("error_ba.txt")),
            stdout           = str(main_output_folder / str("output_ba.txt")),
            inputs           = inputs, 
            outputs          = [bondangle_plot_file], 
            labels           = labels,  
            calib            = calib)

    
    
    #One-shift RP - only for Cd substitutions
    os_data_folder = data_folder / "oneshiftRP"
    os_output_folder = main_output_folder / "oneshiftRP"
    os_output_folder.mkdir(exist_ok=True)

    #oneshfit - Only RP
    os_rp_data_folder = os_data_folder / "onlyRP"
    os_rp_output_folder = os_output_folder / "onlyRP"
    os_rp_output_folder.mkdir(exist_ok=True)
    print("submit one-shift RP simualtions - only RP")
    os_rp_outputtraj_dct = {}
    plane_lst = ["001", "110a", "110b"]
    subst_lst = ["Cs", "Pb4co", "Pb5co", "Pb6co"]
    for plane in plane_lst:
        os_rp_outputtraj_dct[plane] = {}

        #Pure structure
        os_rp_outputtraj_dct[plane]["pure"]={}
        path_atoms = File(str(os_rp_data_folder / str("atoms_"+plane+"_pure.xyz")))
        for sd in range(5):
            simname = plane+"_pure_"+str(sd)
            os_rp_outputtraj_dct[plane]["pure"][sd] = app_run_MD_general(os_rp_output_folder, simname, path_atoms, mult_steps = 2, seed = sd)
            ave_at_apps_lst.append(app_create_ave_struc_general(os_rp_output_folder, os_rp_outputtraj_dct[plane]["pure"][sd], simname))
            if create_input_struc and sd == 0:
                for ind in [0, -1]:
                    new_input_traj.append(app_writeatoms_general(os_rp_output_folder, os_rp_outputtraj_dct[plane]["pure"][sd], simname, index = ind))
        
        #Doped structure
        for subst in subst_lst:
            os_rp_outputtraj_dct[plane][subst] = {}
            for ns in [1,3,6]:
                os_rp_outputtraj_dct[plane][subst][ns] = {}
                path_atoms = File(str(os_rp_data_folder / str("atoms_"+plane+"_Cd_"+subst+"_"+str(ns)+".xyz")))
                for sd in range(5):
                    simname = plane+"_Cd_"+subst+"_"+str(ns)+"_"+str(sd)
                    os_rp_outputtraj_dct[plane][subst][ns][sd] = app_run_MD_general(os_rp_output_folder, simname, path_atoms, mult_steps = 2, seed = sd)
                    ave_at_apps_lst.append(app_create_ave_struc_general(os_rp_output_folder, os_rp_outputtraj_dct[plane][subst][ns][sd], simname))
                    if create_input_struc and sd == 0:
                        for ind in [0, -1]:
                            new_input_traj.append(app_writeatoms_general(os_rp_output_folder,os_rp_outputtraj_dct[plane][subst][ns][sd],simname,index=ind))
        
    #oneshfit -interface
    os_if_data_folder = os_data_folder / "interface"
    os_if_output_folder = os_output_folder / "interface"
    os_if_output_folder.mkdir(exist_ok=True)
    print("submit one-shift RP simualtions - interface")
    #7 different input struc - do MD:
    sim_lst = ["001RP_Cd", "110RP_verysmall", "110RP_verysmall_Cd8", "110RP_small", "110RP_small_Cd8","110RP_normal", "110RP_normal_Cd8"]
    forinput_lst = ["001RP_Cd", "110RP_verysmall", "110RP_verysmall_Cd8", "110RP_small", "110RP_small_Cd8"]
    os_if_outputtraj_dct = {}
    for simname in sim_lst:
        path_atoms = File(str(os_if_data_folder / str("atoms_"+simname+".xyz")))
        os_if_outputtraj_dct[simname] = app_run_MD_general(os_if_output_folder, simname, path_atoms, mult_steps = 2)
        ave_at_apps_lst.append(app_create_ave_struc_general(os_if_output_folder, os_if_outputtraj_dct[simname], simname))
        if create_input_struc and simname in forinput_lst:
            for ind in [0, -1]:
                new_input_traj.append(app_writeatoms_general(os_if_output_folder, os_if_outputtraj_dct[simname], simname, index = ind))


    #big and varying structures, that will need a restart
    for nm in ["_Cd32", ""]:
        simname = "110RP_big_varying" + nm
        path_atoms = File(str(os_if_data_folder / str("atoms_110RP_big_varying"+ nm +".xyz")))
        os_if_outputtraj_dct[simname] = app_run_MD_general(os_if_output_folder, simname, path_atoms)
        ave_at_apps_lst.append(app_create_ave_struc_general(os_if_output_folder, os_if_outputtraj_dct[simname], simname))
        if restart:
            for nr in range(num_res):
                path_atoms = os_if_outputtraj_dct[simname]
                simname = "110RP_big_varying" + nm + "_res" + str(nr)
                os_if_outputtraj_dct[simname] = app_run_MD_general(os_if_output_folder, simname, path_atoms)
                ave_at_apps_lst.append(app_create_ave_struc_general(os_if_output_folder, os_if_outputtraj_dct[simname], simname))

    #Get Pb-I-Pb bond angles distribution
    bondangle_plot_file_1s = File(str(main_output_folder / str("Bondangles_comparison_oneshift.pdf")))
    labels = ["pure_perovskite", "RP_110_2shift", "RP_110_1shift", "interface_110_2shift", "interface_110_1shift"]
    if restart:
        inputs = [puregamma_pv_outputtraj, 
                rp_outputtraj_dct["Pb"]["110"]["Pb"]["Pb6co"][0],
                os_rp_outputtraj_dct["110b"]["pure"][0],
                if_outputtraj_dct["110RP_big_varying_res" + str(num_res - 1)],
                os_if_outputtraj_dct["110RP_big_varying_res" + str(num_res - 1)]]   
    else:
        inputs = [puregamma_pv_outputtraj, 
                rp_outputtraj_dct["Pb"]["110"]["Pb"]["Pb6co"][0],
                os_rp_outputtraj_dct["110b"]["pure"][0],
                if_outputtraj_dct["110RP_big_varying"],
                os_if_outputtraj_dct["110RP_big_varying"]]
    if os.path.exists(bondangle_plot_file_1s.filepath) == False:
        ba_app_1s = app_get_bondangles_distribution(
            execution_folder = main_output_folder,
            stderr           = str(main_output_folder / str("error_ba1s.txt")),
            stdout           = str(main_output_folder / str("output_ba1s.txt")),
            inputs           = inputs, 
            outputs          = [bondangle_plot_file_1s], 
            labels           = labels,  
            calib            = calib)
    
    if create_input_struc:
        #Make trajectory file with all input struc
        path_input_traj = File(str(main_output_folder / str("input_traj.xyz")))
        if os.path.exists(path_input_traj.filepath) == False:
            wit_run = True
            wit_app = app_writeinputtraj(
                execution_folder = main_output_folder,
                stderr           = str(main_output_folder / str("error_wid.txt")),
                stdout           = str(main_output_folder / str("output_wid.txt")),
                inputs           = new_input_traj, 
                outputs          = [path_input_traj])
    
    print("get results of apps")
    parsl.wait_for_current_tasks()

