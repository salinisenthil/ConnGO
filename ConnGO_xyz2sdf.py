import os, sys
import linecache
import numpy as np


input_sdf   = sys.argv[1]
input_xyz   = sys.argv[2]
thresh      = sys.argv[3]
output_sdf  = sys.argv[4]
output_com  = sys.argv[5]
rmsd_name   = sys.argv[6]
charge      =  sys.argv[7]
multiplicity =  sys.argv[8]
mem         = sys.argv[9]
nproc       = sys.argv[10]
#def convert_sdf(input_sdf, input_xyz, thresh, output_sdf, output_com, rmsd_name):

with open(input_sdf, "r") as i_sdf:
        num_lines_isdf = sum(1 for line_isdf in i_sdf)
        P = linecache.getline(input_sdf, 4)
        P_lsp = P.split()
        num_atoms_sdf = int(P_lsp[0])
        count_ext_sdf = int(P_lsp[1])
        R1=np.zeros((num_atoms_sdf,3))
        R3=np.zeros((num_atoms_sdf,3))
        conn_arr=np.zeros((count_ext_sdf,3),dtype=int)
        dist1 = np.zeros(count_ext_sdf)
        dist2 = np.zeros(count_ext_sdf)
        Sy = []

        with open(output_sdf, "w") as new_sdf_f:
            for new_sdf_l in range(1,5):
                pl = linecache.getline(input_sdf, new_sdf_l)
                new_sdf_f.write(pl)

        S1 = []
        for l_isdf in range(5, num_atoms_sdf+5):
            l1 = linecache.getline(input_sdf, l_isdf)
            lsp1 = l1.split()
            ster_l = l1.split(" ",14)
            S1.append(ster_l[14].strip())
            cart1=np.zeros(3)
            cart1 = np.asarray([ float(lsp1[0]) , float(lsp1[1]), float(lsp1[2]) ], dtype=float)
            R1[l_isdf-5][0:3]=cart1
            Sy.append(lsp1[3])

        i = 0
        C1= []
        C2= []
        #for l_isdf1 in range(5+ num_atoms_sdf, 5+num_atoms_sdf +count_ext_sdf):
        for l_isdf1 in range(5+ num_atoms_sdf, num_lines_isdf-1):
            l2 = linecache.getline(input_sdf, l_isdf1)
            C1.append(l2.split())
            C2.append(l2.strip())
            #C1.append(l2.strip())
            lsp2 =  l2.split()
            conn1 = [ int(lsp2[0]), int(lsp2[1]) , int(lsp2[2])  ]
            conn_arr[i][0:3]=conn1
            #print(i,conn1)
            R12=np.zeros(3)
            R12=R1[conn1[0]-1][0:3] - R1[conn1[1]-1][0:3]
            dR12=R12[0]**2 + R12[1]**2 + R12[2]**2
            dist1[i]= round(np.sqrt(dR12),4)
            i = i + 1

with open(output_com,"w") as new_com_f:
    new_com_f.write("%mem="+mem+"gb\n%nproc="+nproc"+\n")
    new_com_f.write("#P METHOD Geom=Connectivity\n\n")
    new_com_f.write(output_com + "\n\n")
    new_com_f.write(charge + "  "+ multiplicity +"\n")
    with open(rmsd_name, "w") as new_rmsd_file:
        with open(input_xyz, "r") as i_xyz:
            num_lines_ixyz = sum(1 for line_ixyz in i_xyz)
            num_atoms_ixyz = linecache.getline(input_xyz, 1)

            with open(output_sdf,"a") as new_sdf_f:
                for l_i in range(3, 3+int(num_atoms_ixyz)):
                    li3 = linecache.getline(input_xyz, l_i)
                    new_com_f.write(li3)
                    lisp3 = li3.split()
                    sdf_newline1 = "    "+ lisp3[1] + "    " +lisp3[2] + "    " + lisp3[3] + " " + Sy[l_i-3] + "   " + S1[l_i-3] + "\n"
                    new_sdf_f.write(sdf_newline1)
                    cart2 = np.zeros(3)
                    cart2 = np.asarray( [ float(lisp3[1]) , float(lisp3[2]), float(lisp3[3]) ], dtype=float  )
                    R3[l_i-3][0:3] = cart2
                for tmp in range(len(C1)):
                    #new_sdf_f.write(C1[tmp] + "\n")
                    new_sdf_f.write(C2[tmp] + "\n")
                new_sdf_f.write("M  END\n$$$$")
                #print(C1[1][0], C1[1][3])
                full =[]
                
                s =1
                for q in range(len(C1)):
                    for k in range(len(C1[q])):
                        full.append(C1[q][k])
                for q in range(len(full)):
                    full[q] = int(full[q])
                print(max(full))



                new_com_f.write("\n")
                #for no in range(len(C1)):
                for no in range(max(full)):
                    E=[]
                    E.append(no+1)
                    #for tmp2 in range(max(full)):
                    for tmp2 in range(len(C1)):
                        if  int(C1[tmp2][0]) == no+1:
                            E.append(int(C1[tmp2][1]))
                            #E.append(int(C1[tmp2][3]))
                            #E.append(float(C1[tmp2][6]))
                            E.append(float(C1[tmp2][2]))
                    E_str = "  "
                    for e in range(len(E)):
                        E_str = E_str + "  "+ str(E[e])
                    #print(E_str)
                    new_com_f.write(E_str + "\n")
                new_com_f.write("\n\n\n")


        nl1=1
        nl2=-2
        for i in range(0, count_ext_sdf):
            R13 = np.zeros(3)
            R13=R3[conn_arr[i][0]-1][0:3] - R3[conn_arr[i][1]-1][0:3]
            dR13=R13[0]**2 + R13[1]**2 + R13[2]**2
            dist2[i]= round(np.sqrt(dR13),4)

            str1='s'
            str2='s'
            if ( dist1[i] > float(thresh)):
                str1='l'
                nl1=nl1*0  
            if ( dist2[i] > float(thresh)):
                str2='l'
                nl2=nl2*0  

            if ( conn_arr[i][2] == 1):
                #print(conn_arr[i][0],conn_arr[i][1],conn_arr[i][2],Sy[conn_arr[i][0]-1],"-",Sy[conn_arr[i][1]-1], dist1[i], str1, dist2[i], str2)
                rm1 = str(conn_arr[i][0]) +"  "+ str(conn_arr[i][1]) +"  "+str(conn_arr[i][2])  +"  "+ str(Sy[conn_arr[i][0]-1])  +"  -  "+ str(Sy[conn_arr[i][1]-1]) +"  "+ str(dist1[i]) +"  "+ str(str1) +"  "+ str(dist2[i]) +"  "+ str(str2) + "\n"
                new_rmsd_file.write(rm1)
            elif ( conn_arr[i][2] == 2):
                rm2 =str(conn_arr[i][0]) +"  "+ str(conn_arr[i][1]) +"  "+ str(conn_arr[i][2]) +"  "+ str(Sy[conn_arr[i][0]-1]) +"  =  "+ str(Sy[conn_arr[i][1]-1]) +"  "+  str(dist1[i]) +"  "+  str(str1) +"  "+  str(dist2[i]) +"  "+  str(str2) + "\n"
                new_rmsd_file.write(rm2)
            elif ( conn_arr[i][2] == 3):
                rm3 = str(conn_arr[i][0]) +"  "+ str(conn_arr[i][1]) +"  "+ str(conn_arr[i][2]) +"  "+ str(Sy[conn_arr[i][0]-1]) +"  #  "+ str(Sy[conn_arr[i][1]-1]) +"  "+ str(dist1[i]) +"  "+  str(str1) +"  "+  str(dist2[i]) +"  "+  str(str2) + "\n"
                new_rmsd_file.write(rm3)
            else:
                new_rmsd_file.write("ERROR: Unknown bond order encountered in SDF!")

#        print(dist2)
        f1_check = 0
        if ( nl1+nl2 == -1):
            new_rmsd_file.write('Both files are OK\n')
        elif ( nl1+nl2 == -2):
            new_rmsd_file.write('File-1 contains long bond!\n')
            f1_check = 1
        elif ( nl1+nl2 == 1):
            new_rmsd_file.write('File-2 contains long bond!\n')
        elif ( nl1+nl2 == 0):
            new_rmsd_file.write('Both files contain long bond!\n')
        RMSD=np.sqrt(sum((dist1-dist2)**2)/count_ext_sdf)
        MaxAD=np.max(np.abs(dist1-dist2))
        MPAD=100*sum(abs((dist1-dist2)/dist1))/count_ext_sdf
        new_rmsd_file.write("RMSD= " + str(round(RMSD,4)) + "\n")
        new_rmsd_file.write("MaxAD= " + str(round(MaxAD,4)) + "\n")
        if "tier1_vs_tier2" in rmsd_name:   # For tier2 - the threshold criteria is MPAD ONLY. MaxAD is waived. if you want to change this, modify this if statement
            if MPAD < 5.0:
                new_rmsd_file.write("MPAD= " +str(round(MPAD,4)) + "  PASS  \n")
            else:
                if f1_check == 1:
                    new_rmsd_file.write("MPAD= " +str(round(MPAD,4)) + "  PASS  \n")
                else:
                    new_rmsd_file.write("MPAD= " +str(round(MPAD,4)) + "  FAIL  \n")
        else:  # For all other tiers, joint threshold of MaxAD and MPAD is enforced.
            if MaxAD < 0.2 and MPAD  < 5.0:
                new_rmsd_file.write("MPAD= " +str(round(MPAD,4)) + "  PASS  \n")
            else:
                if f1_check == 1:
                    new_rmsd_file.write("MPAD= " +str(round(MPAD,4)) + "  PASS  \n")
                else:
                    new_rmsd_file.write("MPAD= " +str(round(MPAD,4)) + "  FAIL  \n")
        


