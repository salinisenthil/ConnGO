import os, sys, string
import linecache, math
import numpy as np


current_dir_path = os.getcwd()
user_file = sys.argv[1]

control_file = current_dir_path + "/control.inp"
job_file = current_dir_path + "/job.inp"


settings = ["load_SDF", "SDF_filename", "Forcefield", "Opt_setting", "Steps", "Convergence_setting", "tier2_mem", "tier2_nprocs","tier2_Functional", "tier2_Basis","tier2_gaussian_setting", "tier3_mem", "tier3_nprocs", "tier3_Functional", "tier3_Basis","tier3_gaussian_setting", "tier4_mem", "tier4_nprocs","tier4_Functional", "tier4_Basis","tier4_gaussian_setting", "charge", "multiplicity", "long_bond", "job_type"]
### READING CONTROL FILE
value = {}
with open(control_file, "r") as ctrl:
    u = 0
    while u < len((settings)):
        for line in ctrl:
            if settings[u] in line:
                val_line = line.split(":")
                if "#" in val_line[1]:
                    val_split = val_line[1].split("#")
                    value[settings[u]] = val_split[0].strip()
                else:
                    value[settings[u]] = val_line[1].strip()
                break
        u = u +1

if value["load_SDF"] == "false":
    initial_smi_path = current_dir_path + "/" + user_file
    smi_fsplit = user_file.split(".")
    smi_file = user_file
    obmin_sdf = smi_fsplit[0] + "_tier1.sdf"
    obmin_xyz = smi_fsplit[0] + "_tier1.xyz"
    
else:
    initial_sdf_path = current_dir_path + "/" + user_file
    obmin_sdf = user_file
    sdf_fsplit = user_file.split(".")
    obmin_xyz = sdf_fsplit[0] + ".xyz"
    smi_file = sdf_fsplit[0] + ".xyz"

def count_pass(filen):
    with open(filen,"r") as cf:
        blc = cf.read()
        if "PASS" in blc:
           count = 1
        else:
           count = 0
    return(count)

if sys.argv[1] == "analyse":
   with open("title_list", "r") as tl:
      total_mol = sum(1 for t_line in tl)
   with open("title_list", "r") as tl:
      tier2_p_count = 0
      tier3_p_count = 0
      tier4_p_count = 0
      for le in tl:
          os.chdir(le.strip())
          t4_check = 0
          if os.path.exists("compare_tier1_vs_tier2"):
              tier2_p_count = tier2_p_count + count_pass("compare_tier1_vs_tier2")
          if os.path.exists("compare_tier1_vs_tier3"):
              tier3_p_count = tier3_p_count + count_pass("compare_tier1_vs_tier3")
          if os.path.exists("compare_tier1_vs_tier4"):
              tier4_p_count = tier4_p_count + count_pass("compare_tier1_vs_tier4")
              t4_check = 1
          if os.path.exists("compare_tier2_vs_tier4"):
              tier4_p_count = tier4_p_count + count_pass("compare_tier2_vs_tier4")
              t4_check = 1
          if os.path.exists("compare_tier3_vs_tier4"):
              tier4_p_count = tier4_p_count + count_pass("compare_tier3_vs_tier4")
              t4_check = 1
          if not (os.path.exists("compare_tier1_vs_tier4") or  os.path.exists("compare_tier2_vs_tier4") or os.path.exists("compare_tier3_vs_tier4")):
              with open("tier4.log","r") as logc:
                  t4out = logc.read()
                  if "Error termination" in t4out  or "imaginary freq" in t4out or t4_check == 0:
                        temp11 = current_dir_path + "/tempfails"
                        with open(temp11, "a") as tmf:
                            tmf.write(le.strip())
                            tmf.write("\n")
                        for files in os.listdir(current_dir_path):
                            if files.endswith(".smi"):
                                smi_file = files
                                com1 = "python " + current_dir_path + "/ConnGO_check_smiles.py  " + le.strip() + "  " + smi_file + "  " + current_dir_path
                                os.system(com1)
                                break
          os.chdir("..")
      if os.path.exists("temp11"):
          with open("temp11", "r") as tempf:
              count_ions = sum(1 for ee in tempf)
      else:
             count_ions = 0
      with open("Results.dat", "w") as res:
          res.write("Total number of molecules entering ConnGO = " + str(total_mol) + "\n\n")
          res.write("TIER-2\n")
          res.write("Number of molecules PASS = " + str(tier2_p_count) + "\n")
          res.write("Number of molecules FAIL = " + str(total_mol-tier2_p_count) + "\n\n")
          res.write("TIER-3\n")
          res.write("Number of molecules PASS = " + str(tier3_p_count) + "\n")
          res.write("Number of molecules FAIL = " + str(total_mol-tier2_p_count-tier3_p_count) + "\n\n")
          res.write("TIER-4\n")
          res.write("Number of molecules PASS = " + str(tier4_p_count) + "\n")
          res.write("Number of molecules FAIL = " + str(total_mol-tier4_p_count) + "\n\n")  ## WRONG
          res.write("Number of zwitterions encountered that fail ConnGO = " + str(count_ions) + "  (if greater than 0, see modified_smiles.smi)\n\n")
          res.write("List of molecules failing ConnGO:\n")
          n = 1
          if os.path.exists("tempfails"):
              for fail_l in open("tempfails", "r"):
                  res.write(str(n) + "  " + fail_l )
                  n = n+1

          os.system("rm -f temp*  ")
      sys.exit()


if sys.argv[1] == "launch":
        if value["job_type"] == "submission":
            #os.system("ls run*.sh > runlist")
            with open("title_list", "r") as tl:
                for linet in tl:
                    os.chdir(linet.strip())
                    pt1 = "launching " + linet.strip()
                    print(pt1)
                    os.system("qsub run.sh")
                    os.chdir("..")
 
        if value["job_type"] == "interactive":
            with open("title_list", "r") as tl:
                for linet in tl:
                    os.chdir(linet.strip())
                    pt1 = "launching " + linet.strip()
                    print(pt1)
                    os.system("bash run.sh  &")
                    os.chdir("..")

        sys.exit()

if value["load_SDF"] == "false":
    #obmin
    to_genxyz = "obabel  -oxyz  " + smi_file + " >  " + smi_fsplit[0] +  "_initial.xyz  --gen3d"
    #to_obmin = "obminimize -osdf -ff MMFF94 -sd -n 100000 -c 1e-8  " + smi_fsplit[0] + "_initial.xyz  >  " + OMST_sdf
    to_obmin = "obminimize  -osdf  -ff  " + value["Forcefield"] + "  -" + value["Opt_setting"]  + "  -n  " + value["Steps"] + "  -c  " + value["Convergence_setting"] +"  " + smi_fsplit[0] + "_initial.xyz  >  " + obmin_sdf
    omst_xyzgen = "obabel -oxyz " + obmin_sdf + "  > " + obmin_xyz
    
    print("generating initial xyz from smi")
    os.system(to_genxyz)
    print("OMST in progess")
    os.system(to_obmin)
    print("OMST done, sdf created")
    os.system(omst_xyzgen)
    print("OMST xyz created")



def launch_test(title_name, sdf_file,xyz_file, tier_level):
        fl1 = tier_level+"_Functional"
        bl1 = tier_level+"_Basis"
        st1 = tier_level+"_gaussian_setting"
        tt1 = tier_level+"_mem"
        tt2 = tier_level+"_nprocs"
        #current_launch_com = current_title + "_" + tier_level + ".com"
        current_launch_com = tier_level + ".com"
        create_comfile = "python  " + current_dir_path +"/ConnGO_xyz2sdf.py   " + sdf_file + "  " + xyz_file + "  1.7  discard.sdf  " + current_launch_com + "  discard_rmsd   " + value["charge"].strip() +"  "+ value["multiplicity"].strip() +"  "+ value[tt1].strip() +"  " + value[tt2].strip() + " \n\n"
        #method_insert_hf =  "sed -i \"s/METHOD/HF\\\/sto-3G SCF(maxcycles=100) Opt(calcall,cartesian,tight,maxcyc=100) Int(Grid=Ultrafine) freq /\"  " + current_launch_com + " \n\n"
        #method_insert_hf =  "sed -i \"s/METHOD/" + func + "\\\/" + basis1 + "  SCF(maxcycles=100) Opt(calcall,cartesian,tight,maxcyc=100) Int(Grid=Ultrafine) freq /\" " + current_launch_com + " \n\n"
        method_insert_hf =  "sed -i \"s/METHOD/" + value[fl1] + "\\\/" + value[bl1] + "  " + value[st1] +"  /\"  " + current_launch_com + " \n\n"

        launch_gaussian = "g16  " + current_launch_com  + "\n\n"

        run_file.write(create_comfile)
        run_file.write(method_insert_hf)
        run_file.write(launch_gaussian)

def hf_fail(title_name):
        launch_test(current_title,current_sdf, current_xyz , "tier3" )
        check_ET_321log = "if grep  -q \"Error termination\"  " + current_321g_log + " ;  then\n"
        run_file.write(check_ET_321log)
        launch_test(current_title, current_sdf,current_xyz, "tier4")
        check_2dfp_status(title_name, current_sdf, "tier1")
        check_IF_321log = "elif grep  -q \"imaginary freq\"  " + current_321g_log + " ;  then\n"
        run_file.write(check_IF_321log)
        launch_test(current_title, current_sdf,current_xyz, "tier4")
        check_2dfp_status(title_name, current_sdf, "tier1")
        run_file.write("else\n")

        opt_321g_xyz = "for f in " + current_321g_log+ "; do NATOMS=$(grep -m1 NAtom $f | awk \'{print $2}\'); echo $NATOMS; echo ${f%.*}; grep -i -A $(( $NATOMS + 4 )) \'Standard Orientation\' $f|  tail -$NATOMS | awk \'{printf \"%3g\" \"%15.8f\" \"%15.8f\" \"%15.8f\", $2, $4, $5, $6}{printf \"\\n\"}\' ; done >  " +current_321g_xyz+ " \n\n"
        #create_321g_rmsd = "/apps/moldis_utils/sdf2sdf/xyz2sdf_v1.exe  " +  current_sdf +"   " + current_321g_xyz +  "  321g.sdf  1.7  " + " temp1.com   2  > rmsd_mmff_dp_" + current_title + "\n\n"
        create_321g_rmsd = "python  " + current_dir_path + "/ConnGO_xyz2sdf.py  " + current_sdf + "  " + current_321g_xyz + "  1.7 " + current_321g_sdf + " temp.com   compare_tier1_vs_tier3   "  + value["charge"].strip() +"  "+ value["multiplicity"].strip() + " 8  2  \n\n"
        check_FAIL_321g = "if grep  -q \"FAIL\"  compare_tier1_vs_tier3   ;  then\n"

        run_file.write(opt_321g_xyz)
        run_file.write(create_321g_rmsd)
        run_file.write(check_FAIL_321g)

        launch_test(current_title,current_sdf, current_xyz , "tier4")
        check_2dfp_status(title_name, current_sdf, "tier1")

        run_file.write("else\n")
        launch_test(current_title, current_321g_sdf,current_321g_xyz, "tier4")
        check_2dfp_status(title_name,current_321g_sdf, "tier3")

        run_file.write(" \n fi \n fi \n")



def check_2dfp_status(title_name, previous_tier_sdf, tier):

    check_normalterm = "if grep  -q \"Normal termination\"  " + current_2dfp_log + " ;  then\n\n"
    
    opt_2dfp_xyz = "for f in " + current_2dfp_log + "; do NATOMS=$(grep -m1 NAtom $f | awk \'{print $2}\'); echo $NATOMS; echo ${f%.*}; grep -i -A $(( $NATOMS + 4 )) \'Standard Orientation\' $f|  tail -$NATOMS | awk \'{printf \"%3g\" \"%15.8f\" \"%15.8f\" \"%15.8f\", $2, $4, $5, $6}{printf \"\\n\"}\' ; done >  " + current_2dfp_xyz + " \n\n"
   
    #create_2dfp_rmsd = "/apps/moldis_utils/sdf2sdf/xyz2sdf_v1.exe  " +  previous_tier_sdf + "  " +current_2dfp_xyz +  "  2dfp.sdf  1.7  " + "temp1.com   2  > rmsd_"+ tier+"_2dfp_" + current_title + "\n\n"
    create_2dfp_rmsd = "python  " + current_dir_path + "/ConnGO_xyz2sdf.py  " + previous_tier_sdf + "  " + current_2dfp_xyz + "  1.7 " + current_2dfp_sdf + "  temp.com  compare_" + tier + "_vs_tier4    " + value["charge"].strip() +"  "+ value["multiplicity"].strip() +" 8  2 \n\n"
    run_file.write(check_normalterm)
    run_file.write(opt_2dfp_xyz)
    run_file.write(create_2dfp_rmsd)
    copy_files(title_name)
    run_file.write("\nfi\n")

def copy_files(title_name):
    if value["job_type"] == "submission":
        copy_back = "cp *com *log compare* *xyz  *sdf  " + current_dir_path + "/"+title_name+"/ \n"
        run_file.write(copy_back)
        run_file.write("rm -rf $WORKDIR\n\n")
        run_file.write("exit\n\n")
    else:
        run_file.write("exit\n\n")



######## Create seperate sdf files from OMST sdf
ln_sdf = 1
with open(obmin_sdf, "r") as initial_sdf:
    num_lines_sdf = sum(1 for line_par in initial_sdf)
    while ln_sdf <= num_lines_sdf:
        P = linecache.getline(obmin_sdf, ln_sdf+3)
        P_lsp = P.split()
        num_atoms_sdf = int(P_lsp[0])
        count_ext_sdf = int(P_lsp[1])
        next_ln_sdf = ln_sdf + num_atoms_sdf + count_ext_sdf + 6
        P1 = linecache.getline(obmin_sdf, ln_sdf).strip()
        fname_sdf = P1 + "_tier1.sdf"
        with open("title_list", "a") as dl1:
             dl1.write(P1.strip())
             dl1.write("\n")
        os.mkdir(P1)
        os.chdir(P1)
        for line_sdf in range(ln_sdf, next_ln_sdf):
            with open("tier1.sdf","a") as new_f_sdf:
                pline_sdf = linecache.getline(obmin_sdf, line_sdf)
                new_f_sdf.write(pline_sdf)
        os.chdir("..")
        ln_sdf = next_ln_sdf

print("seperate OMST sdfs created")
######## Create seperate xyz files from OMST xyz
ln_mxyz = 1
with open(obmin_xyz, "r") as ini_xyz:
    num_lines_mxyz = sum(1 for line_mxyz in ini_xyz)
    countd = 1
    while ln_mxyz <= num_lines_mxyz:
        num_atoms_mxyz = linecache.getline(obmin_xyz, ln_mxyz)
        next_ln_mxyz = ln_mxyz + int(num_atoms_mxyz) + 2
        with open("title_list", "r") as dr:
           dname = linecache.getline("title_list", countd).strip()
        os.chdir(dname)
        with open("tier1.xyz", "w") as new_f_mxyz:
            new_f_mxyz.write(num_atoms_mxyz)
            new_f_mxyz.write(dname)
            new_f_mxyz.write("\n")
            for ln_mxyz in range(ln_mxyz + 2, next_ln_mxyz):
                ln_mxyz = linecache.getline(obmin_xyz, ln_mxyz)
                new_f_mxyz.write(ln_mxyz)
            os.chdir("..")
        countd = countd + 1
        ln_mxyz = next_ln_mxyz

print("seperate OMST xyzs created")


####### Begining of FLOW 3
for line_title in open("title_list","r"):
    current_title = line_title.strip()
    #run_script =  "run_" + current_title + ".sh"
    run_script =  "run.sh"
    os.chdir(current_title)
    with open(run_script, "a") as run_file:
        if value["job_type"] == "submission":
            with open(job_file, "r") as job_block:
                j = job_block.read()
                run_file.write(j)
                run_file.write("\n")

        current_sdf      = "tier1.sdf"
        current_xyz      = "tier1.xyz"
        current_hf_log   = "tier2.log"
        current_hf_sdf   = "tier2.sdf"
        current_hf_xyz   = "tier2.xyz"
        current_321g_log = "tier3.log"
        current_321g_sdf = "tier3.sdf"
        current_321g_xyz = "tier3.xyz"
        current_2dfp_log = "tier4.log"
        current_2dfp_sdf = "tier4.sdf"
        current_2dfp_xyz = "tier4.xyz"


        if value["job_type"] == "submission":
            #copy_smi = "cp   " + current_dir_path+ "/" + smi_file  +  " . \n"
            new_path = os.getcwd()
            copy_sdf = "cp   " + new_path+ "/" + current_sdf  +  " . \n"
            copy_xyz  = "cp   " + new_path+ "/" + current_xyz  +  " .\n\n"
            copy_ctrl  = "cp   " + current_dir_path+ "/control.inp  .\n\n"
            
            #run_file.write(copy_smi)
            run_file.write(copy_sdf)
            run_file.write(copy_xyz)
            run_file.write(copy_ctrl)
        launch_test(current_title, current_sdf, current_xyz , "tier2")

        check_ET_HFlog = "if grep  -q \"Error termination\"  " + current_hf_log + " ;  then\n"
        run_file.write(check_ET_HFlog)
        hf_fail(current_title)
        check_IF_HFlog = "elif grep  -q \"imaginary freq\"  " + current_hf_log + " ;  then\n"
        run_file.write(check_IF_HFlog)
        hf_fail(current_title)
        run_file.write("\nelse\n")

        opt_hf_xyz = " for f in " + current_hf_log+ "; do NATOMS=$(grep -m1 NAtom $f | awk \'{print $2}\'); echo $NATOMS; echo ${f%.*}; grep -i -A $(( $NATOMS + 4 )) \'Standard Orientation\' $f|  tail -$NATOMS | awk \'{printf \"%3g\" \"%15.8f\" \"%15.8f\" \"%15.8f\", $2, $4, $5, $6}{printf \"\\n\"}\' ; done >  " +current_hf_xyz+ " \n\n"
        #create_hf_rmsd = "/apps/moldis_utils/sdf2sdf/xyz2sdf_v1.exe  " +  current_sdf + "  "+ current_hf_xyz +  "  hf.sdf  1.7  " + "temp1.com   2  > rmsd_mmff_hf_" + current_title + "\n\n"
        create_hf_rmsd = "python " + current_dir_path + "/ConnGO_xyz2sdf.py  " + current_sdf + "  " + current_hf_xyz + "  1.7 " + current_hf_sdf +"  temp.com  compare_tier1_vs_tier2   " + value["charge"].strip() +"  "+ value["multiplicity"].strip() +" 8  2\n\n"
        check_FAIL_hf = "if grep  -q \"FAIL\"  compare_tier1_vs_tier2   ;  then\n"
        
        run_file.write(opt_hf_xyz)
        run_file.write(create_hf_rmsd)
        run_file.write(check_FAIL_hf)
        hf_fail(current_title)
        run_file.write("else\n")

        launch_test(current_title, current_hf_sdf, current_hf_xyz, "tier4")
        check_2dfp_status(current_title, current_hf_sdf, "tier2")
        run_file.write("\n fi \n fi\n")
        print("created "+ run_script)
    os.chdir("..")
