import os, sys, string
import linecache, math
import numpy as np


current_dir_path = os.getcwd()
smi_file = sys.argv[1]

initial_smi_path = current_dir_path + "/" + smi_file
control_file = current_dir_path + "/control.inp"
job_file = current_dir_path + "/job.inp"

smi_fsplit = smi_file.split(".")
obmin_sdf = smi_fsplit[0] + "_obmin.sdf"
obmin_xyz = smi_fsplit[0] + "_obmin.xyz"

settings = ["Forcefield", "Opt_setting", "Steps", "Convergence_setting", "tier2_Functional", "tier2_Basis","tier2_gaussian_setting", "tier3_Functional", "tier3_Basis","tier3_gaussian_setting", "tier4_Functional", "tier4_Basis","tier4_gaussian_setting", "charge", "multiplicity", "long_bond", "job_type"]
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

#print(value["multiplicity"])

def count_occ(t_filename):
    count=0
    for line in open(t_filename,"r"):
        name = line.strip()
        with open(name, "r") as rf:
            blc = rf.read()
            if "PASS" in blc:
                count = count + 1
    return(count)

if sys.argv[1] == "analyse":
    # tier2 vs 4 , 3 vs 4 and 1 vs 4 - 2dfp  # 1 vs 2 - hf    # 1 vs 3 - dp 
    os.system("ls rmsd_tier2_vs_tier4* > atemp1")
    lc1 = sum(1 for ll1 in open("atemp1","r"))
    if lc1 == 0:
        t2_4_pass = 0
    else:
        t2_4_pass = count_occ("atemp1")

    os.system("ls rmsd_tier3_vs_tier4* > atemp2")
    lc2 = sum(1 for ll2 in open("atemp2","r"))
    if lc2 == 0:
        t3_4_pass = 0
    else:
        t3_4_pass = count_occ("atemp2")

    os.system("ls rmsd_tier1_vs_tier4* > atemp3")
    lc3 = sum(1 for ll3 in open("atemp3","r"))
    if lc3 == 0:
        t1_4_pass = 0
    else:
        t1_4_pass = count_occ("atemp3")

    os.system("ls rmsd_tier1_vs_tier2* > atemp4")
    lc4 = sum(1 for ll4 in open("atemp4","r"))
    if lc4 == 0:
        t1_2_pass = 0
    else:
        t1_2_pass = count_occ("atemp4")

    os.system("ls rmsd_tier1_vs_tier3* > atemp5")
    lc5 = sum(1 for ll5 in open("atemp5","r"))
    if lc5 == 0:
        t1_3_pass = 0
    else:
        t1_3_pass = count_occ("atemp5")

    
    with open("title_list", "r") as tl:
        num_mol = sum(1 for t_line in tl)
    if os.path.exists("modified_tautomers.smi"):
        with open("modified_tautomers.smi", "r") as mt:
            num_mod_taut = sum(1 for mt_line in mt)
    else:
        num_mod_taut = 0

    os.system("grep -l FAIL rmsd*tier4*  > atemp6")
    for aline1 in open("atemp6", "r"):
        if "tier2" in aline1:
            m1 = aline1.replace("rmsd_tier2_vs_tier4_","")
        if "tier3" in aline1:
            m1 = aline1.replace("rmsd_tier3_vs_tier4_","")
        if "tier1" in aline1:
            m1 = aline1.replace("rmsd_tier1_vs_tier4_","")
        with open("temp_fails", "a") as tmp_f:
            tmp_f.write(m1.strip() + "\n")
    os.system("grep -l \"Error termination\"  *tier4.log  > atemp7")
    os.system("grep -l \"imaginary freq\"  *tier4.log  >> atemp7")
    for aline2 in open("atemp7", "r"):
        if "log" in aline2:
            m2 = aline2.replace("_tier4.log","")
        with open("temp_fails", "a") as tmp_f:
            tmp_f.write(m2.strip() + "\n")
        

    with open("Results.dat", "w") as res:
        res.write("Total number of molecules entering ConnGO = " + str(num_mol) + "\n\n")
        res.write("TIER-2\n")
        res.write("Number of molecules PASS = " + str(t1_2_pass) + "\n")
        res.write("Number of molecules FAIL = " + str(num_mol-t1_2_pass) + "\n\n")
        res.write("TIER-3\n")
        res.write("Number of molecules PASS = " + str(t1_3_pass) + "\n")
        res.write("Number of molecules FAIL = " + str(num_mol-t1_2_pass-t1_3_pass) + "\n\n")
        res.write("TIER-4\n")
        res.write("Number of molecules PASS = " + str(t2_4_pass+t3_4_pass+t1_4_pass) + "\n")
        res.write("Number of molecules FAIL = " + str(num_mol-(t2_4_pass+t3_4_pass+t1_4_pass)) + "\n\n")  ## WRONG
        res.write("Number of tautomers encountered that fail ConnGO = " + str(num_mod_taut) + "  (if greater than 0, see modified_tautomers.smi)\n\n")
        res.write("List of molecules failing ConnGO:\n")
        n = 1
        for fail_l in open("temp_fails", "r"):
            res.write(str(n) + "  " + fail_l )
            n = n+1
        
    os.system("rm -f atemp* temp* discard* ")
    sys.exit()

if sys.argv[1] == "launch":
        if value["job_type"] == "submission":
            os.system("ls run*.sh > runlist")
            for rl in open("runlist", "r"):
                launch = "qsub " + rl
                pr = "launching " + rl
                print(pr)
                os.system(launch)
 
        if value["job_type"] == "interactive":
            os.system("ls run*.sh > runlist")
            for rl in open("runlist", "r"):
                 launch = "bash " + rl + " & "
                 pr = "launching " + rl
                 print(pr)
                 os.system(launch)
        os.system("rm -f runlist")
        sys.exit()

#obmin
to_genxyz = "obabel  -oxyz  " + smi_file + " >  " + smi_fsplit[0] +  "_initial.xyz  --gen3d"
to_obmin = "obminimize -osdf -ff MMFF94 -sd -n 100000 -c 1e-8  " + smi_fsplit[0] + "_initial.xyz  >  " + OMST_sdf
to_obmin = "obminimize  -osdf  -ff  " + value["Forcefield"] + "  -" + value["Opt_setting"]  + "  -n  " + value["Steps"] + "  -c  " + value["Convergence_setting"] +"  " + smi_fsplit[0] + "_initial.xyz  >  " + OMST_sdf
omst_xyzgen = "obabel -oxyz " + OMST_sdf + "  > " + OMST_xyz

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
        current_launch_com = current_title + "_" + tier_level + ".com"
        #create_comfile = "/apps/moldis_utils/sdf2sdf/xyz2sdf_v1.exe  " + current_sdf + "  " + xyz_file + "   discard.sdf  1.7  " + current_launch_com + " 2  > discard_rmsd" + current_title + "\n\n"
        create_comfile = "python  " + current_dir_path +"/ConnGO_xyz2sdf.py   " + sdf_file + "  " + xyz_file + "  1.7  discard.sdf  " + current_launch_com + "  discard_rmsd   " + value["charge"].strip() +"  "+ value["multiplicity"].strip() +" \n\n"
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
        create_321g_rmsd = "python  " + current_dir_path + "/ConnGO_xyz2sdf.py  " + current_sdf + "  " + current_321g_xyz + "  1.7 " + current_321g_sdf + " temp.com   rmsd_tier1_vs_tier3_" + title_name+ "   " + value["charge"].strip() +"  "+ value["multiplicity"].strip() + " \n\n"
        check_FAIL_321g = "if grep  -q \"FAIL\"  rmsd_tier1_vs_tier3_" + current_title + "  ;  then\n"

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

    check_IF_2dfplog = "if grep  -q \"imaginary freq\"  " + current_2dfp_log + " ;  then\n"
    check_taut = "python " + current_dir_path + "/ConnGO_check_tautomer.py  " + title_name + "  " + smi_file + "  " + current_dir_path + " \n"
    run_file.write(check_IF_2dfplog)
    run_file.write(check_taut)
    copy_files()

    check_normalterm = "elif grep  -q \"Normal termination\"  " + current_2dfp_log + " ;  then\n\n"
    
    opt_2dfp_xyz = "for f in " + current_2dfp_log + "; do NATOMS=$(grep -m1 NAtom $f | awk \'{print $2}\'); echo $NATOMS; echo ${f%.*}; grep -i -A $(( $NATOMS + 4 )) \'Standard Orientation\' $f|  tail -$NATOMS | awk \'{printf \"%3g\" \"%15.8f\" \"%15.8f\" \"%15.8f\", $2, $4, $5, $6}{printf \"\\n\"}\' ; done >  " + current_2dfp_xyz + " \n\n"
   
    #create_2dfp_rmsd = "/apps/moldis_utils/sdf2sdf/xyz2sdf_v1.exe  " +  previous_tier_sdf + "  " +current_2dfp_xyz +  "  2dfp.sdf  1.7  " + "temp1.com   2  > rmsd_"+ tier+"_2dfp_" + current_title + "\n\n"
    create_2dfp_rmsd = "python  " + current_dir_path + "/ConnGO_xyz2sdf.py  " + previous_tier_sdf + "  " + current_2dfp_xyz + "  1.7 " + current_2dfp_sdf + "  temp.com  rmsd_" + tier + "_vs_tier4_" + title_name + "   " + value["charge"].strip() +"  "+ value["multiplicity"].strip() +"\n\n"
    run_file.write(check_normalterm)
    run_file.write(opt_2dfp_xyz)
    run_file.write(create_2dfp_rmsd)
    check_FAIL_2dfp = "if grep  -q \"FAIL\"  rmsd_"+tier+"_vs_tier4_" + title_name + " ;  then\n"
    check_taut1 = "python  "+current_dir_path+"/ConnGO_check_tautomer.py  " + title_name + "  " + smi_file + "  " + current_dir_path + "\n"
    run_file.write(check_FAIL_2dfp)
    run_file.write(check_taut1)
    run_file.write("\nfi\n")
    copy_files()


    run_file.write("\nelse\n")
    check_taut = "python  " + current_dir_path + "/ConnGO_check_tautomer.py  " + title_name + "  " + smi_file + "  " + current_dir_path + "\n"
    run_file.write(check_taut)
    copy_files()
    run_file.write("\nfi\n")

def copy_files():
    if value["job_type"] == "submission":
        copy_back = "cp *com *log rmsd* *xyz  *sdf  " + current_dir_path + "/\n"
        run_file.write(copy_back)
        run_file.write("rm -rf $WORKDIR\n\n")
        run_file.write("exit\n\n")
    else:
        run_file.write("exit\n\n")



######## Create seperate sdf files from OMST sdf
os.system("rm *tier1.sdf")
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
        for line_sdf in range(ln_sdf, next_ln_sdf):
            with open(fname_sdf,"a") as new_f_sdf:
                pline_sdf = linecache.getline(obmin_sdf, line_sdf)
                new_f_sdf.write(pline_sdf)
        ln_sdf = next_ln_sdf

print("seperate OMST sdfs created")

######## Create seperate xyz files from OMST xyz
ln_mxyz = 1
with open(obmin_xyz, "r") as ini_xyz:
    num_lines_mxyz = sum(1 for line_mxyz in ini_xyz)
    while ln_mxyz <= num_lines_mxyz:
        num_atoms_mxyz = linecache.getline(obmin_xyz, ln_mxyz)
        next_ln_mxyz = ln_mxyz + int(num_atoms_mxyz) + 2
        fln_mxyz = linecache.getline(obmin_xyz, ln_mxyz + 1).strip()
        flsp_mxyz = fln_mxyz.split()
        f_name_mxyz = flsp_mxyz[0] + "_tier1.xyz"
        with open(f_name_mxyz, "w") as new_f_mxyz:
            new_f_mxyz.write(num_atoms_mxyz)
            new_f_mxyz.write(f_name_mxyz)
            new_f_mxyz.write("\n")
            for ln_mxyz in range(ln_mxyz + 2, next_ln_mxyz):
                ln_mxyz = linecache.getline(obmin_xyz, ln_mxyz)
                new_f_mxyz.write(ln_mxyz)
        ln_mxyz = next_ln_mxyz

print("seperate OMST xyzs created")

os.system("ls *tier1.xyz > title_list")
os.system("sed -i 's/_tier1.xyz//g' title_list")


####### Begining of FLOW 3
for line_title in open("title_list","r"):
    current_title = line_title.strip()
    run_script =  "run_" + current_title + ".sh"
    with open(run_script, "a") as run_file:
        if value["job_type"] == "submission":
            with open(job_file, "r") as job_block:
                j = job_block.read()
                run_file.write(j)
                run_file.write("\n")

        current_sdf = current_title + "_tier1.sdf"
        current_xyz = current_title + "_tier1.xyz"
        current_hf_log = current_title + "_tier2.log"
        current_hf_sdf = current_title + "_tier2.sdf"
        current_hf_xyz = current_title + "_tier2.xyz"
        current_321g_log = current_title + "_tier3.log"
        current_321g_sdf = current_title + "_tier3.sdf"
        current_321g_xyz = current_title + "_tier3.xyz"
        current_2dfp_log = current_title + "_tier4.log"
        current_2dfp_sdf = current_title + "_tier4.sdf"
        current_2dfp_xyz = current_title + "_tier4.xyz"

        if value["job_type"] == "submission":
            #copy_smi = "cp   " + current_dir_path+ "/" + smi_file  +  " . \n"
            copy_sdf = "cp   " + current_dir_path+ "/" + current_sdf  +  " . \n"
            copy_xyz  = "cp   " + current_dir_path+ "/" + current_xyz  +  " .\n\n"
            
            #run_file.write(copy_smi)
            run_file.write(copy_sdf)
            run_file.write(copy_xyz)
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
        create_hf_rmsd = "python " + current_dir_path + "/ConnGO_xyz2sdf.py  " + current_sdf + "  " + current_hf_xyz + "  1.7 " + current_hf_sdf +"  temp.com  rmsd_tier1_vs_tier2_" + current_title + "  " + value["charge"].strip() +"  "+ value["multiplicity"].strip() +" \n\n"
        check_FAIL_hf = "if grep  -q \"FAIL\"  rmsd_tier1_vs_tier2_" + current_title + "  ;  then\n"
        
        run_file.write(opt_hf_xyz)
        run_file.write(create_hf_rmsd)
        run_file.write(check_FAIL_hf)
        hf_fail(current_title)
        run_file.write("else\n")

        launch_test(current_title, current_hf_sdf, current_hf_xyz, "tier4")
        check_2dfp_status(current_title, current_hf_sdf, "tier2")
        run_file.write("\n fi \n fi\n")
        print("created "+ run_script)
