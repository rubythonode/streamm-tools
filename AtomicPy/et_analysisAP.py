#! /usr/bin/env python
# Read in group ref file and pull data from qm calc

# Dr. Travis Kemper
# NREL
# Initial Date 12/18/2013
# Email travis.kemper@nrel.gov
# Version 2.00 
#


def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] "
    parser = OptionParser(usage=usage)

    parser.add_option("-v","--verbose", dest="verbose", default=False, help="Verbose output ")

    parser.add_option("--in_et", dest="in_et", type="string", default="et.dat", help="Iinput file of V_AB values  ")

    # Gromacs
    parser.add_option("--in_top", dest="in_top", type="string", default="", help="Input gromacs topology file ")
    parser.add_option("--in_gro", dest="in_gro", type="string", default="", help="Input gromacs structure file ")
 
    # Groups
    parser.add_option("--group_ptma", dest="group_ptma", default=True, help="Group TEMPO molecules ")
    parser.add_option("--nlist_bonds", dest="nlist_bonds", default=True, help="Build neighbor list from bonds")
    parser.add_option("--atom_types", dest="atom_types", type="string", default="", help="Read atom types that will replace default elements ")

    (options, args) = parser.parse_args()

    return options, args

def ket(vij):
    import math 

    dg=0
    lamb=1.02
    kb=8.6173324*10**-5
    hbar=6.58211928*10**-16
    
    T=298.0

    ket = vij**2/hbar*math.sqrt( math.pi/(kb*T*lamb))*math.e**( -(lamb+dg)**2/(4*lamb*kb*T))
    
    return ket

def sigma_m(N,ave,ave_sq):
    """
    Calculate the standard deviation of the mean for a confidence
    interval of 90%. Will return zero for single valued data sets 
    """
    import numpy
    
    # Some website that probably does not exist
    #   http://mathworld.wolfram.com/Studentst-Distribution.html
    #   http://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm
    #
    
    v = N - 1 #  Degrees of freedom
    
    # Set appropriate Students t prefactor 
    #if( v > 100 ):
    #	Coefint_pre = 1.64487
    #el
    
    if( v > 30 ):
	Coefint_pre = 1.66023
    elif( v > 10 ):
	Coefint_pre = 1.69726
    elif( v > 5 ):
	Coefint_pre = 1.81246
    elif( v == 5 ):
	Coefint_pre = 2.01505
    elif( v == 4 ):
	Coefint_pre = 2.13185
    elif( v == 3 ):
	Coefint_pre = 2.35336
    elif( v == 2 ):
	Coefint_pre = 2.91999
    elif( v == 1 ):
	Coefint_pre = 6.31375  
	
    if( N > 1 ):
	v_sqrt = numpy.sqrt(  N - 1 )
	sigma = numpy.sqrt(  ( ave_sq ) - (ave)**2 ) # Standard deviation
	sigma_m = Coefint_pre*sigma/v_sqrt
    else:
	sigma_m = 0.0  # Set zero for unknow error 
    
    return sigma_m

def main():
    import os, sys,  math
    import numpy  as np 
    import file_io, gromacs , top, elements, groups, prop
    
    debug = 0

    FOURPI =  4.0*math.pi
    VOL_CONST = FOURPI/3.0

    E_cg = -480.099900537626/0.03674931
    E_ng = -480.347576682435/0.03674931

    print "E_ng ",E_ng
    print "E_cg ",E_cg
 
    # Get options
    options, args = get_options()
    
    # Get reference structure file 
    if( len(options.in_gro) ):
        if( options.verbose ): print "  Reading in ",options.in_gro
        GTYPE, R, VEL, LV = gromacs.read_gro(options.in_gro)
        
    # Read in topology file
    if( len(options.in_top) ):
        if( options.verbose ): print "  Reading in ",options.in_top
        ATYPE , RESN , RESID , GTYPE ,CHARN , CHARGES ,AMASS,BONDS,ANGLES,DIH, MOLNUMB, MOLPNT, MOLLIST  = gromacs.read_top(options.in_top)
        ASYMB , ELN  = elements.mass_asymb(AMASS) 
    
    # Retype special atom types and replace element
    if( len( options.atom_types)): 
        if( options.verbose ): print "  Reading in ",options.atom_types
        ASYMB , ELN  = top.special_types(ATYPE,ASYMB , ELN , options.atom_types)
        
    # Print system information
    if( options.verbose ):
	print " prop "
        #prop.print_prop( AMASS,ELN,LV,CHARGES )
        ##top.print_prop( BONDS,ANGLES,DIH )
        
    # Create neighbor list
    if( options.nlist_bonds ):
        if( options.verbose ): print "  Creating neighbor list with bonds"
        NBLIST,NBINDEX  = groups.build_nablist_bonds(ELN,BONDS)

    # Find groups
    if( options.group_ptma ):
        if( options.verbose ): print "  Find groups of TEMPO "
        group_index_i,group_list_i,group_numb,group_cnt = groups.tempo(  ATYPE,ELN,NBLIST,NBINDEX, options )

    # Read in values for pairs of groups 
    et_file = options.in_et

    if( file_io.file_exists(et_file)):
        f_et_in = open(et_file,"r")
        Lines = f_et_in.readlines()
        f_et_in.close()

    pair_cnt = 0 
    site_cnt = 0

    for line in Lines:
        col = line.split()
        if ( len(col) > 1 and col[0].strip() != "#" ):
            pair_cnt += 1 
            site_i = int( col[0] )
            if( site_i > site_cnt ): site_cnt = site_i
            site_i = int( col[1])
            if( site_i > site_cnt ): site_cnt = site_i
    
    # Initialize histogram settings

    en_r_list = np.zeros( site_cnt )

    en_c_list = np.zeros( site_cnt )
    en_rc_list = np.zeros( site_cnt )

    I_rp_list = np.zeros( pair_cnt )
    I_pr_list = np.zeros( pair_cnt )

    r_bin_size = 0.1 #(angstroms)
    r_min = 0.0
    r_max = 10.0
    r_floor = int(r_min/r_bin_size)*r_bin_size 
    r_range = r_max - r_floor
    n_r_bins = int(r_range/r_bin_size) + 1
    vol_cut = VOL_CONST*r_max*r_max*r_max

    r_bin_cnt = np.zeros(n_r_bins)

    et_bins = np.zeros(n_r_bins)
    et_sq_bins  = np.zeros(n_r_bins)
    et_qu_bins  = np.zeros(n_r_bins)

    et_sq_bins_intra  = np.zeros(n_r_bins)
    et_qu_bins_intra  = np.zeros(n_r_bins)
    r_bin_cnt_intra = np.zeros(n_r_bins)
    et_sq_bins_inter  = np.zeros(n_r_bins)
    et_qu_bins_inter  = np.zeros(n_r_bins)
    r_bin_cnt_inter = np.zeros(n_r_bins)

    dG_r_bins  = np.zeros(n_r_bins)
    dGsq_r_bins  = np.zeros(n_r_bins)
    I_r_bins  = np.zeros(n_r_bins)
    Isq_r_bins  = np.zeros(n_r_bins)
    dG_r_bins_intra  = np.zeros(n_r_bins)
    dGsq_r_bins_intra  = np.zeros(n_r_bins)
    I_r_bins_intra  = np.zeros(n_r_bins)
    Isq_r_bins_intra  = np.zeros(n_r_bins)
    dG_r_bins_inter  = np.zeros(n_r_bins)
    dGsq_r_bins_inter  = np.zeros(n_r_bins)
    I_r_bins_inter  = np.zeros(n_r_bins)
    Isq_r_bins_inter  = np.zeros(n_r_bins)

    lambda_r_bins  = np.zeros(n_r_bins)
    lambdasq_r_bins  = np.zeros(n_r_bins)
    lambda_r_bins_intra  = np.zeros(n_r_bins)
    lambdasq_r_bins_intra  = np.zeros(n_r_bins)
    lambda_r_bins_inter  = np.zeros(n_r_bins)
    lambdasq_r_bins_inter  = np.zeros(n_r_bins)

    dG_bin_size = 0.2
    dG_min = -1.0 
    dG_max = 1.0

    dG_floor = int(dG_min/dG_bin_size)*dG_bin_size 
    dG_range = dG_max - dG_floor
    dG_n_bins = int(dG_range/dG_bin_size) + 1

    print "dG_n_bins",dG_n_bins

    dG_r_dG_bins  = np.zeros((n_r_bins,dG_n_bins))
    dG_r_dG_cnt  = np.zeros((n_r_bins,dG_n_bins))
    dG_r_dG_bins_intra  = np.zeros((n_r_bins,dG_n_bins))
    dG_r_dG_cnt_intra  = np.zeros((n_r_bins,dG_n_bins))
    dG_r_dG_bins_inter  = np.zeros((n_r_bins,dG_n_bins))
    dG_r_dG_cnt_inter  = np.zeros((n_r_bins,dG_n_bins))

    Vsq_sum = 0.0
    Vsq_rsq_sum = 0.0
    Vsq_rqu_sum = 0.0
              
    Vsq_sum_inter = 0.0
    Vsq_rsq_sum_inter = 0.0
                
    Vsq_sum_intra = 0.0
    Vsq_rsq_sum_intra = 0.0
                        
    newdat = open("newdat.dat","w")
    
    grp_et = np.empty((3000,3000),dtype="float")
    
    print len(grp_et)
    print len(grp_et[0])
    print grp_et[0][0]

    # Find volume
    volume_i =   prop.volume( LV ) 

    print "volume_i",volume_i

    pair_indx = -1 

    for line in Lines:
        col = line.split()
        if( col[0] != "#" ):
            pair_indx += 1

            n_i = int( col[0] )   # Nuetral group #
            c_i = int( col[1] )   # Cation group #
            dr =  float(col[3] )  # inter nitrogen spacing


            E_i_r = float( col[4] ) 
            E_i_p = float( col[5] ) 
            E_j_r = float( col[6] ) 
            E_j_p = float( col[7] )

            E_ij_rp = float( col[8] )
            E_ij_pr = float( col[9] )
            et_value  = float( col[10] )

            lambda_inner_ij = E_i_r - E_ng + E_j_p - E_cg
            lambda_inner_ji = E_j_r - E_ng + E_i_p - E_cg

            delta_Gij = E_ij_pr - E_ij_rp
            delta_Gji =   E_ij_rp - E_ij_pr

            I_rp = E_ij_rp - E_i_r - E_j_p
            I_pr = E_ij_pr - E_i_p - E_j_r

            I_rp_list[pair_indx] = abs(I_rp)
            I_pr_list[pair_indx] = abs(I_pr)
        
            r_bin_v = int( dr/ r_bin_size )   # Bin value 
            if( r_bin_v < n_r_bins ):         
                r_bin_cnt[r_bin_v] += 1
		#

                check_O = True 
                
                # Find molecule numbers 
                Ni_o = group_index_i[n_i]
                atom_i = group_list_i[Ni_o]
                mol_i = MOLNUMB[atom_i]
                r_i = np.array( R[atom_i] )

                if( check_O ):

                    N_o = group_index_i[n_i]
                    N_f = group_index_i[n_i + 1 ] - 1
                    atom_k = -1

                    # print " checking mol ",n_i,N_o,N_f,group_list_i[N_o],group_list_i[N_f]
                    # print " PTMA group i ",n_i
                    
                    for a_indx in range( N_o,N_f):
                        test_k = group_list_i[a_indx]
                        if( ATYPE[test_k] == "ON" ):
                            atom_k = test_k
                    if( atom_k > 0 ):
                        # print " atom_k ",atom_k,ATYPE[atom_k],mol_i,MOLNUMB[atom_k]
                        # print "index ",atom_i," or index ",atom_k
                        r_k = np.array( R[atom_k] )
                    else:
                        error_line = " No Oxygen found in molecule %d "%(mol_i)
                        sys.exit(error_line)

                # Find molecule numbers 
                Nj_o = group_index_i[c_i]
                atom_j = group_list_i[Nj_o]
                mol_j = MOLNUMB[atom_j]
                r_j = np.array( R[atom_j] )

                if( check_O ):

                    N_o = group_index_i[c_i]
                    N_f = group_index_i[c_i + 1 ] - 1
                    atom_l = -1

                    # print " checking mol ",n_i,N_o,N_f,group_list_i[N_o],group_list_i[N_f]
                    # print " PTMA group j ",c_i
                    
                    for a_indx in range( N_o,N_f):
                        test_k = group_list_i[a_indx]
                        if( ATYPE[test_k] == "ON" ):
                            atom_l = test_k
                    if( atom_l > 0 ):
                        # print " atom_k ",atom_k,ATYPE[atom_k],mol_i,MOLNUMB[atom_k]
                        # print "index ",atom_i," or index ",atom_l
                        r_l = np.array( R[atom_l] )
                    else:
                        error_line = " No Oxygen found in molecule %d "%(mol_i)
                        sys.exit(error_line)

                r_ij =  prop.mag_drij_c(r_i,r_j,LV) 

                r_error = abs(r_ij - dr )
                
                if( r_error > 1.0 ):
                    sys.exit(" bad r_ij value ")


                #
                # r_i --- r_k
                #  | theta_kij
                #  |
                #  | theta_ijl 
                #  r_j --- r_l
                #
                #   r_k   phi_kijl    r_l
                #        \         /
                #           r_ij
                #

                r_ik = prop.rij_c(r_i,r_k,LV)
                r_ij = prop.rij_c(r_i,r_j,LV)
                r_ji = prop.rij_c(r_j,r_i,LV)
                r_jl = prop.rij_c(r_j,r_l,LV)
                
                theta_kij = prop.getAngle(r_ik,r_ij)
                theta_ijl = prop.getAngle(r_ji,r_jl)
                dih_lijl = prop.getDihedral(r_k,r_i,r_j,r_l,LV)

                debug_theta = False
                if( debug_theta ):

                    print "index ",atom_i," or index ",atom_k," or index ",atom_j," or index ",atom_l

                    print " r_i ",r_i
                    print " r_j ",r_j
                    print " r_k ",r_k
                    print " r_l ",r_l
                    print " r_ik ",r_ik
                    print " r_ij ",r_ij
                    print " r_ji ",r_ji
                    print " r_jl ",r_jl
                    print " theta_kij ",theta_kij
                    print " theta_ijl ",theta_ijl
                    print " dih_lijl ",dih_lijl

                    sys.exit(" theta ")

                newdat.write(" %d %d %f %f %f %f %f  %f %f  %f %f  %f %f \n"% (  n_i , c_i , dr , E_i_r , E_i_p, E_j_r, E_j_p, E_ij_rp, E_ij_pr, et_value, theta_kij, theta_ijl, dih_lijl))

                    
                et_bins[r_bin_v] += et_value
                et_sq = et_value*et_value 
                et_qu = et_sq*et_sq 
                et_sq_bins[r_bin_v] += et_sq
                et_qu_bins[r_bin_v] += et_qu
                grp_et[n_i][c_i] = et_value

                dG_r_bins[r_bin_v] += abs( delta_Gij )
                dGsq_r_bins[r_bin_v] += abs( delta_Gij*delta_Gij )
                I_r_bins[r_bin_v] += abs(I_rp)
                Isq_r_bins[r_bin_v] += abs(I_rp*I_rp)
                I_r_bins[r_bin_v] += abs(I_pr)
                Isq_r_bins[r_bin_v] += abs(I_pr*I_pr)

                lambda_r_bins[r_bin_v] += abs( lambda_inner_ij )
                lambdasq_r_bins[r_bin_v] += abs( lambda_inner_ij*lambda_inner_ij )
                lambda_r_bins[r_bin_v] += abs( lambda_inner_ji )
                lambdasq_r_bins[r_bin_v] += abs( lambda_inner_ji*lambda_inner_ji )

                dGij_bin_v = int( ( delta_Gij - dG_floor ) / dG_bin_size )   # Bin value
                dG_r_dG_bins[r_bin_v,dGij_bin_v] += 1
                dGji_bin_v = int( ( delta_Gji - dG_floor ) / dG_bin_size )   # Bin value
                dG_r_dG_bins[r_bin_v,dGji_bin_v] += 1


                Vsq_sum += et_sq
                Vsq_rsq_sum += et_sq*dr*dr
                Vsq_rqu_sum +=  et_sq*dr*dr*dr*dr

                if(  mol_i == mol_j):
                    et_sq_bins_intra[r_bin_v] += et_sq
                    et_qu_bins_intra[r_bin_v] += et_qu
                    r_bin_cnt_intra[r_bin_v] += 1


                    dG_r_bins_intra[r_bin_v] += abs( delta_Gij )
                    I_r_bins_intra[r_bin_v] += abs(I_rp)
                    I_r_bins_intra[r_bin_v] += abs(I_pr)

                    dGsq_r_bins_intra[r_bin_v] += abs( delta_Gij*delta_Gij )
                    Isq_r_bins_intra[r_bin_v] += abs(I_rp*I_rp)
                    Isq_r_bins_intra[r_bin_v] += abs(I_pr*I_pr)


                    lambda_r_bins_intra[r_bin_v] += abs( lambda_inner_ij )
                    lambdasq_r_bins_intra[r_bin_v] += abs( lambda_inner_ij*lambda_inner_ij )
                    lambda_r_bins_intra[r_bin_v] += abs( lambda_inner_ji )
                    lambdasq_r_bins_intra[r_bin_v] += abs( lambda_inner_ji*lambda_inner_ji )

                    dGij_bin_v = int( ( delta_Gij - dG_floor ) / dG_bin_size )   # Bin value
                    dG_r_dG_bins_intra[r_bin_v,dGij_bin_v] += 1
                    dGji_bin_v = int( ( delta_Gji - dG_floor ) / dG_bin_size )   # Bin value
                    dG_r_dG_bins_intra[r_bin_v,dGji_bin_v] += 1

                    Vsq_sum_intra += et_sq
                    Vsq_rsq_sum_intra += et_sq*dr*dr
                    
                else:
                    et_sq_bins_inter[r_bin_v] += et_sq
                    et_qu_bins_inter[r_bin_v] += et_qu
                    r_bin_cnt_inter[r_bin_v] += 1

                    dG_r_bins_inter[r_bin_v] += abs( delta_Gij )
                    I_r_bins_inter[r_bin_v] += abs(I_rp)
                    I_r_bins_inter[r_bin_v] += abs(I_pr)

                    dGsq_r_bins_inter[r_bin_v] += abs( delta_Gij*delta_Gij )
                    Isq_r_bins_inter[r_bin_v] += abs(I_rp*I_rp)
                    Isq_r_bins_inter[r_bin_v] += abs(I_pr*I_pr)

                    lambda_r_bins_inter[r_bin_v] += abs( lambda_inner_ij )
                    lambdasq_r_bins_inter[r_bin_v] += abs( lambda_inner_ij*lambda_inner_ij )
                    lambda_r_bins_inter[r_bin_v] += abs( lambda_inner_ji )
                    lambdasq_r_bins_inter[r_bin_v] += abs( lambda_inner_ji*lambda_inner_ji )

                    dGij_bin_v = int( ( delta_Gij - dG_floor ) / dG_bin_size )   # Bin value
                    dG_r_dG_bins_inter[r_bin_v,dGij_bin_v] += 1
                    dGji_bin_v = int( ( delta_Gji - dG_floor ) / dG_bin_size )   # Bin value
                    dG_r_dG_bins_inter[r_bin_v,dGji_bin_v] += 1

                    
                    Vsq_sum_inter += et_sq
                    Vsq_rsq_sum_inter += et_sq*dr*dr

		if( debug ):
		    if( r_bin_v == 50 ):
			print dr , et_value

    newdat.close()
    
    print " group_cnt ",group_cnt
    print " pair_indx ",pair_indx
    total_pairs = pair_indx + 1
    box_vol_ave = volume_i


    l_sq = Vsq_rsq_sum/Vsq_sum
    print " r^2 rhogv / rhogv ",l_sq," l ",np.sqrt(l_sq)

    l_sq = Vsq_rqu_sum/Vsq_rsq_sum
    print " sum r^4 rhogv / sum r^2 rhogv ",l_sq," l ",np.sqrt(l_sq)

    l_sq_intra = Vsq_rsq_sum_intra/Vsq_sum_intra
    print " l^2 intra ",l_sq_intra," l ",np.sqrt(l_sq_intra)

    l_sq_inter = Vsq_rsq_sum_inter/Vsq_sum_inter
    print " l^2 inter ",l_sq_inter," l ",np.sqrt(l_sq_inter)
                       
    et_cut = 0.0000001
    for g_i in range( 500):
        for g_j in range( 500):
            if( grp_et[g_i][g_j] > et_cut and   grp_et[g_j][g_i] > et_cut  ):
                print g_i,g_j,grp_et[g_i][g_j] ,  grp_et[g_j][g_i] 
         
    #
    # Calaculate averages, g_v, 4 pi 
    #
    #  b - bin number
    #  N - number of groups
    # dV - volume element  (A^3)
    #
    #  < V_{AB} > (b)   = sum(VAB)(b) / cnt(b) 
    #  < V_{AB}^2 > (b) = sum(VAB^2)(b) / cnt(b) 
    #
    #             dVg_v = 2 * sum(VAB^2) / N 
    #            dr_vol = dV 
    #            rhog_v = dVg_v/dr_vol            # (rho * g_v )
    #   
    #  
    
    
    Coef_interval =  1.00 # 96

    SUM_vsq = 0.0 
    SUM_rsq_vsq = 0.0 
    SUM_rqu_vsq = 0.0 
    
    SUM_rhog_v = 0.0
    SUM_rsq_rhog_v = 0.0
    SUM_rsq_rhog_v_sig_mean = 0.0
    SUM_rqd_rhog_v = 0.0
    SUM_rqd_rhog_v_sig_mean = 0.0
    SUM_rhog_v_intra = 0.0
    SUM_rsq_rhog_v_intra = 0.0
    SUM_rsq_rhog_v_intra_sig_mean = 0.0
    SUM_rqd_rhog_v_intra = 0.0
    SUM_rqd_rhog_v_intra_sig_mean = 0.0
    SUM_rhog_v_inter = 0.0
    SUM_rsq_rhog_v_inter = 0.0
    SUM_rsq_rhog_v_inter_sig_mean = 0.0
    SUM_rqd_rhog_v_inter = 0.0
    SUM_rqd_rhog_v_inter_sig_mean = 0.0
    
    f_gr = open("gr.hist","w")
    
    f_VAB = open("VAB.hist","w")
    f_VABsq = open("VABsq.hist","w")
    f_rquVABsq = open("rquVABsq.hist","w")
    f_dG_r = open("dG_r.hist","w")
    f_I_r = open("I_r.hist","w")
    f_lambda_r = open("Lambda_inner_r.hist","w")

    f_dG_r3D = open("dG_r3D.hist","w")

    
    #f_gv = open("rhogv.hist","w")
    #f_rsqgv = open("rsqgv.hist","w")
    #f_rqugv = open("rqugv.hist","w")
    #f_dVgv = open("dVgv.hist","w")

    f_gr.write("# r_ij ; cnts ; g(r) ")

    f_VAB.write("# r_ij ; cnts ; <V_{AB}> (eV);  mstd( <V_{AB}>(eV)  ) , intra  < V_{AB} >  (eV); intra mstd( <V_{AB}> (eV)  ), inter  < V_{AB} >  (eV); inter mstd( <V_{AB}>(eV)   ")
    f_VABsq.write("# r_ij ; cnts ; < V_{AB}^2 >  (eV^2); mstd( <V_{AB}^2>(eV^2)  ),   intra  < V_{AB}^2 >  (eV^2); intra mstd( <V_{AB}^2>(eV^2)  ), inter  < V_{AB}^2 >  (eV^2); inter mstd( <V_{AB}^2>(eV^2)  ) ")
    f_rquVABsq.write("# r_ij ; cnts ; < V_{AB}^2 >  (eV^2); mstd( <V_{AB}^2>(eV^2)  ),   intra  < V_{AB}^2 >  (eV^2); intra mstd( <V_{AB}^2>(eV^2)  ), inter  < V_{AB}^2 >  (eV^2); inter mstd( <V_{AB}^2>(eV^2)  ) ")

    f_dG_r.write("# r_ij ; cnts ; <dG> (eV);  mstd( <dG>(eV)  ) , intra  <dG>  (eV); intra mstd( <dG> (eV)  ), inter  <dG >  (eV); inter mstd( <dG>(eV)   ")
    f_I_r.write("# r_ij ; cnts ; <I> (eV);  mstd( <I>(eV)  ) , intra  <I>  (eV); intra mstd( <I> (eV)  ), inter  <I >  (eV); inter mstd( <I>(eV)   ")
    f_lambda_r.write("# r_ij ; cnts ; <lambda> (eV);  mstd( <lambda>(eV)  ) , intra  <lambda>  (eV); intra mstd( <lambda> (eV)  ), inter  <lambda>  (eV); inter mstd( <lambda>(eV)   ")
    
    # f_gv.write("# r_ij ; cnts ; <rhog_v> (eV^2/A^3);  mstd( <rhog_v>(eV^2/A^3)  ); intra <rhog_v> (eV^2/A^3); intra  mstd( <rhog_v>(eV^2/A^3)  ); inter <rhog_v> (eV^2/A^3); inter mstd( <rhog_v>(eV^2/A^3) )  ")
    # f_et_h = open("et_h.hist","w")
    # f_et.write("# r_ij ; cnts ;  sum V_{AB} ; < V_{AB} >  ; sum V_{AB}^2 ; < V_{AB}^2 > ; g_v ; < V_{AB}^2 >_intra ; g_v_intra ; < V_{AB}^2 >_inter ; g_v_inter ; 4 pi r^2 g_v")
    # f_VAB.write("# r_ij ; cnts ; <V_{AB}>(eV);  mstd( <V_{AB}>(eV)  )  ")
    # f_VABsq.write("# r_ij ; cnts ; < V_{AB}^2 >  (eV^2); mstd( <V_{AB}^2>(eV^2)  ),   intra  < V_{AB}^2 >  (eV^2); intra mstd( <V_{AB}^2>(eV^2)  ), inter  < V_{AB}^2 >  (eV^2); inter mstd( <V_{AB}^2>(eV^2)  ) ")
    # f_rsqgv.write("# r_ij ; cnts ; <rhog_v> (eV^2/A);  mstd( <rhog_v>(eV^2/A)  ); intra <rhog_v> (eV^2/A); intra  mstd( <rhog_v>(eV^2/A)  ); inter <rhog_v> (eV^2/A); inter mstd( <rhog_v>(eV^2/A) )  ")
    # f_rqugv.write("# r_ij ; cnts ; <rhog_v> (eV^2 A);  mstd( <rhog_v>(eV^2 A)  ); intra <rhog_v> (eV^2 A); intra  mstd( <rhog_v>(eV^2 A)  ); inter <rhog_v> (eV^2 A); inter mstd( <rhog_v>(eV^2 A) )  ")
    # f_dVgv.write("# r_ij ; cnts ; <dVrhog_v> (eV^2);  mstd( <dVrhog_v> (eV^2)  ); intra <dVrhog_v> (eV^2); intra  mstd(<dVrhog_v> (eV^2)); inter <dVrhog_v> (eV^2); inter mstd( <dVrhog_v> (eV^2) )  ")
    # r_ij,r_bin_cnt[b_indx],b_sum, b_average, et_sq ,et_sq_ave, rhog_v ,et_sq_intra_ave,rhog_v_intra ,et_sq_inter_ave, rhog_v_inter

    total_pairs = site_cnt*(site_cnt -1 ) # 
    total_cnts = np.sum( r_bin_cnt )

    print "total_cnts",total_cnts

    inv_rho_box_ij = float(box_vol_ave)/float(total_pairs)

    n_shperes  = site_cnt
    inv_rho_sphere_ij = vol_cut*float(n_shperes)/float(total_cnts)

    print "inv_rho_box_ij ",inv_rho_box_ij,box_vol_ave,total_pairs
    print "inv_rho_sphere_ij ",inv_rho_sphere_ij,vol_cut,float(n_shperes),float(total_cnts)

    for b_indx in range(n_r_bins):
        r_ij = b_indx*r_bin_size + r_floor + r_bin_size*0.5
        r_in = r_ij - r_bin_size*0.5
        r_out = r_ij + r_bin_size*0.5
        dr_vol = VOL_CONST*( r_out**3 - r_in**3 )

        N_cnts = r_bin_cnt[b_indx]

        # g(r)
        box_gr  = 0.0 # float(N_cnts)/dr_vol 

        # Average Vij
        et_sum = 0.0 #  et_bins[b_indx] 
        et_average = 0.0 # et_sum/float(N_cnts)
        et_sq = 0.0 # et_sq_bins[b_indx]
        et_sq_ave = 0.0 # et_sq/float(N_cnts)
        # Average Vij^2
        et_qu = 0.0 # et_qu_bins[b_indx]
        et_qu_ave = 0.0 # et_qu/float(N_cnts)

        et_sig_mean = 0.0 # sigma_m(N_cnts,et_average,et_sq_ave)
        et_sq_sig_mean = 0.0 # sigma_m(N_cnts,et_sq_ave,et_qu_ave)

        rqu_et_sq_ave = 0.0 #  dr_qu*et_sq_ave
        rqu_et_sq_ave_sig_mean = 0.0 # dr_qu*et_sq_sig_mean

        dG_ave = 0.0 # dG_r_bins[b_indx]/float(N_cnts)
        dGsq_ave = 0.0 # dGsq_r_bins[b_indx]/float(N_cnts)
        dG_ave_sig_mean = 0.0 # sigma_m(N_cnts,dG_ave,dGsq_ave)


        I_ave = 0.0 # I_r_bins[b_indx]/float(N_cnts*2)
        Isq_ave = 0.0 # Isq_r_bins[b_indx]/float(N_cnts*2)
        I_ave_sig_mean = 0.0 # sigma_m(N_cnts*2,I_ave,Isq_ave)


        lambda_ave = 0.0 # lambda_r_bins[b_indx]/float(N_cnts*2)
        lambdasq_ave = 0.0 # lambdasq_r_bins[b_indx]/float(N_cnts*2)
        lambda_ave_sig_mean = 0.0 # sigma_m(N_cnts*2,lambda_ave,lambdasq_ave)

        et_sq_intra = 0.0 # et_sq_bins_intra[b_indx]
        et_sq_intra_ave = 0.0 # et_sq_intra/N_cnt_intra
        et_qu_intra = 0.0 # et_qu_bins_intra[b_indx]
        et_qu_intra_ave = 0.0 # et_qu_intra/N_cnt_intra

        et_sq_sig_mean_intra = 0.0 # sigma_m(N_cnt_intra,et_sq_intra_ave,et_qu_intra_ave)

        rqu_et_sq_ave_intra = 0.0 #  dr_qu*et_sq_intra_ave
        rqu_et_sq_ave_sig_mean_intra = 0.0 # dr_qu*et_sq_sig_mean_intra

        dG_ave_intra = 0.0 # dG_r_bins_intra[b_indx]/float(N_cnts)
        dGsq_ave_intra = 0.0 # dGsq_r_bins_intra[b_indx]/float(N_cnts)
        dG_ave_sig_mean_intra = 0.0 # sigma_m(N_cnts,dG_ave_intra,dGsq_ave_intra)

        I_ave_intra = 0.0 # I_r_bins_intra[b_indx]/float(N_cnts*2)
        Isq_ave_intra = 0.0 # Isq_r_bins_intra[b_indx]/float(N_cnts*2)
        I_ave_sig_mean_intra = 0.0 # sigma_m(N_cnts*2,I_ave_intra,Isq_ave_intra)



        lambda_ave_intra = 0.0 # lambda_r_bins_intra[b_indx]/float(N_cnts*2)
        lambdasq_ave_intra = 0.0 # lambdasq_r_bins_intra[b_indx]/float(N_cnts*2)
        lambda_ave_sig_mean_intra = 0.0 # sigma_m(N_cnts*2,lambda_ave_intra,lambdasq_ave_intra)


        et_sq_inter = 0.0 # et_sq_bins_inter[b_indx]
        et_sq_inter_ave = 0.0 # et_sq_inter/N_cnt_inter
        et_qu_inter = 0.0 # et_qu_bins_inter[b_indx]
        et_qu_inter_ave = 0.0 # et_qu_inter/N_cnt_inter

        et_sq_sig_mean_inter = 0.0 # sigma_m(N_cnt_inter,et_sq_inter_ave,et_qu_inter_ave)


        rqu_et_sq_ave_inter = 0.0 #  dr_qu*et_sq_inter_ave
        rqu_et_sq_ave_sig_mean_inter = 0.0 # dr_qu*et_sq_sig_mean_inter

        dG_ave_inter = 0.0 # dG_r_bins_inter[b_indx]/float(N_cnts)
        dGsq_ave_inter = 0.0 # dGsq_r_bins_inter[b_indx]/float(N_cnts)
        dG_ave_sig_mean_inter = 0.0 # sigma_m(N_cnts,dG_ave_inter,dGsq_ave_inter)

        I_ave_inter = 0.0 # I_r_bins_inter[b_indx]/float(N_cnts*2)
        Isq_ave_inter = 0.0 # Isq_r_bins_inter[b_indx]/float(N_cnts*2)
        I_ave_sig_mean_inter = 0.0 # sigma_m(N_cnts*2,I_ave_inter,Isq_ave_inter)

        lambda_ave_inter = 0.0 # lambda_r_bins_inter[b_indx]/float(N_cnts*2)
        lambdasq_ave_inter = 0.0 # lambdasq_r_bins_inter[b_indx]/float(N_cnts*2)
        lambda_ave_sig_mean_inter = 0.0 # sigma_m(N_cnts*2,lambda_ave_inter,lambdasq_ave_inter)

        if(  N_cnts > 0 ):
	    
	    dr_sq = r_ij*r_ij
	    dr_qu = dr_sq*dr_sq

            # g(r)
            box_gr  = float(N_cnts)/dr_vol 

	    # Average Vij
            et_sum =  et_bins[b_indx] 
            et_average = et_sum/float(N_cnts)
            et_sq = et_sq_bins[b_indx]
            et_sq_ave = et_sq/float(N_cnts)
            # Average Vij^2
            et_qu = et_qu_bins[b_indx]
            et_qu_ave = et_qu/float(N_cnts)

	    et_sig_mean = sigma_m(N_cnts,et_average,et_sq_ave)
	    et_sq_sig_mean = sigma_m(N_cnts,et_sq_ave,et_qu_ave)

            rqu_et_sq_ave =  dr_qu*et_sq_ave
            rqu_et_sq_ave_sig_mean = dr_qu*et_sq_sig_mean
            
            dG_ave = dG_r_bins[b_indx]/float(N_cnts)
            dGsq_ave = dGsq_r_bins[b_indx]/float(N_cnts)
            dG_ave_sig_mean = sigma_m(N_cnts,dG_ave,dGsq_ave)

            
            I_ave = I_r_bins[b_indx]/float(N_cnts*2)
            Isq_ave = Isq_r_bins[b_indx]/float(N_cnts*2)
            I_ave_sig_mean = sigma_m(N_cnts*2,I_ave,Isq_ave)


            lambda_ave = lambda_r_bins[b_indx]/float(N_cnts*2)
            lambdasq_ave = lambdasq_r_bins[b_indx]/float(N_cnts*2)
            lambda_ave_sig_mean = sigma_m(N_cnts*2,lambda_ave,lambdasq_ave)

            #g_v = et_sq/dr_vol/float( 2424 )

            SUM_vsq +=  et_sq 
            SUM_rsq_vsq +=  dr_sq*et_sq 
            SUM_rqu_vsq +=  dr_qu*et_sq 

            # Old crap
	    
            # rho g_v(r)
            # dVg_v =  2.0*et_sq/float( group_cnt )
            # rhog_v = dVg_v/dr_vol
	    # dVg_v_sig_mean = 2.0*et_sq_sig_mean*float(N_cnts)/float( group_cnt )	    
	    # rhog_v_sig_mean = dVg_v_sig_mean/dr_vol

            
	    # fourpi_rsqr_rhog_v = FOURPI*dr_sq*rhog_v
	    # fourpi_rsqr_rhog_v_sig_mean = FOURPI*dr_sq*rhog_v_sig_mean 

            # rqu_rhog_v = dr_qu*rhog_v
	    # rqu_rhog_v_sig_mean =dr_qu*rhog_v_sig_mean 

            # SUM_rhog_v += rhog_v
            # SUM_rsq_rhog_v += dr_sq*rhog_v
	    # SUM_rsq_rhog_v_sig_mean  += dr_sq*rhog_v_sig_mean*dr_sq*rhog_v_sig_mean
            # SUM_rqd_rhog_v += dr_qu*rhog_v
	    # SUM_rqd_rhog_v_sig_mean += dr_qu*rhog_v_sig_mean*dr_qu*rhog_v_sig_mean

	    N_cnt_intra = float( r_bin_cnt_intra[b_indx] )
	    N_cnt_inter = float( r_bin_cnt_inter[b_indx] )
            
            if(  N_cnt_intra > 0 ):
                    
                et_sq_intra = et_sq_bins_intra[b_indx]
                et_sq_intra_ave = et_sq_intra/N_cnt_intra
                et_qu_intra = et_qu_bins_intra[b_indx]
                et_qu_intra_ave = et_qu_intra/N_cnt_intra
		
		et_sq_sig_mean_intra = sigma_m(N_cnt_intra,et_sq_intra_ave,et_qu_intra_ave)

                rqu_et_sq_ave_intra =  dr_qu*et_sq_intra_ave
                rqu_et_sq_ave_sig_mean_intra = dr_qu*et_sq_sig_mean_intra



                dG_ave_intra = dG_r_bins_intra[b_indx]/float(N_cnt_intra)
                dGsq_ave_intra = dGsq_r_bins_intra[b_indx]/float(N_cnt_intra)
                dG_ave_sig_mean_intra = sigma_m(N_cnt_intra,dG_ave_intra,dGsq_ave_intra)


                I_ave_intra = I_r_bins_intra[b_indx]/float(N_cnt_intra*2)
                Isq_ave_intra = Isq_r_bins_intra[b_indx]/float(N_cnt_intra*2)
                I_ave_sig_mean_intra = sigma_m(N_cnt_intra*2,I_ave_intra,Isq_ave_intra)



                lambda_ave_intra = lambda_r_bins_intra[b_indx]/float(N_cnt_intra*2)
                lambdasq_ave_intra = lambdasq_r_bins_intra[b_indx]/float(N_cnt_intra*2)
                lambda_ave_sig_mean_intra = sigma_m(N_cnt_intra*2,lambda_ave_intra,lambdasq_ave_intra)


            
                #dVg_v_intra = 2.0*et_sq_intra/float( group_cnt )
		#dVg_v_sig_mean_intra = 2.0*et_sq_sig_mean_intra*N_cnt_intra/float( group_cnt )
                #rhog_v_intra = dVg_v_intra/dr_vol
		#rhog_v_sig_mean_intra = dVg_v_sig_mean_intra/dr_vol
		
		
		#fourpi_rsqr_rhog_v_intra = FOURPI*dr_sq*rhog_v_intra
		#fourpi_rsqr_rhog_v_sig_mean_intra = FOURPI*dr_sq*rhog_v_sig_mean_intra 
	    
                #SUM_rhog_v_intra += rhog_v_intra
                #SUM_rsq_rhog_v_intra += dr_sq*rhog_v_intra
		#SUM_rsq_rhog_v_intra_sig_mean += dr_sq*rhog_v_sig_mean_intra*dr_sq*rhog_v_sig_mean_intra
                #SUM_rqd_rhog_v_intra += dr_qu*rhog_v_intra
                #SUM_rqd_rhog_v_intra_sig_mean += dr_qu*rhog_v_sig_mean_intra*dr_qu*rhog_v_sig_mean_intra
                
            if(  N_cnt_inter > 0 ):
                et_sq_inter = et_sq_bins_inter[b_indx]
                et_sq_inter_ave = et_sq_inter/N_cnt_inter
                et_qu_inter = et_qu_bins_inter[b_indx]
                et_qu_inter_ave = et_qu_inter/N_cnt_inter
		
		et_sq_sig_mean_inter = sigma_m(N_cnt_inter,et_sq_inter_ave,et_qu_inter_ave)

		
                rqu_et_sq_ave_inter =  dr_qu*et_sq_inter_ave
                rqu_et_sq_ave_sig_mean_inter = dr_qu*et_sq_sig_mean_inter

                dG_ave_inter = dG_r_bins_inter[b_indx]/float(N_cnt_inter)
                dGsq_ave_inter = dGsq_r_bins_inter[b_indx]/float(N_cnt_inter)
                dG_ave_sig_mean_inter = sigma_m(N_cnt_inter,dG_ave_inter,dGsq_ave_inter)

                I_ave_inter = I_r_bins_inter[b_indx]/float(N_cnt_inter*2)
                Isq_ave_inter = Isq_r_bins_inter[b_indx]/float(N_cnt_inter*2)
                I_ave_sig_mean_inter = sigma_m(N_cnt_inter*2,I_ave_inter,Isq_ave_inter)

                lambda_ave_inter = lambda_r_bins_inter[b_indx]/float(N_cnt_inter*2)
                lambdasq_ave_inter = lambdasq_r_bins_inter[b_indx]/float(N_cnt_inter*2)
                lambda_ave_sig_mean_inter = sigma_m(N_cnt_inter*2,lambda_ave_inter,lambdasq_ave_inter)


                # dVg_v_inter = 2.0*et_sq_inter/float( group_cnt )
		# dVg_v_sig_mean_inter = 2.0*et_sq_sig_mean_inter*N_cnt_inter/float( group_cnt )
                # rhog_v_inter = dVg_v_inter/dr_vol
		# rhog_v_sig_mean_inter = dVg_v_sig_mean_inter/dr_vol
		
		# fourpi_rsqr_rhog_v_inter = FOURPI*dr_sq*rhog_v_inter
		# fourpi_rsqr_rhog_v_sig_mean_inter = FOURPI*dr_sq*rhog_v_sig_mean_inter 
	    
		
                # SUM_rhog_v_inter += rhog_v_inter
                # SUM_rsq_rhog_v_inter += dr_sq*rhog_v_inter
                # SUM_rsq_rhog_v_inter_sig_mean += dr_sq*rhog_v_sig_mean_inter*dr_sq*rhog_v_sig_mean_inter
                # SUM_rqd_rhog_v_inter += dr_qu*rhog_v_inter
                # SUM_rqd_rhog_v_inter_sig_mean += dr_qu*rhog_v_sig_mean_inter*dr_qu*rhog_v_sig_mean_inter
                
            
            
        f_gr.write( "\n %f %d %16.12f "% (r_ij , r_bin_cnt[b_indx], box_gr ) )
        f_VABsq.write( "\n %f %d %16.12f %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , r_bin_cnt[b_indx], et_sq_ave, et_sq_sig_mean, et_sq_intra_ave, et_sq_sig_mean_intra, et_sq_inter_ave, et_sq_sig_mean_inter  ) )
        f_rquVABsq.write( "\n %f %d  %16.12f %16.12f %16.12f  %16.12f  %16.12f  %16.12f "% (r_ij , r_bin_cnt[b_indx], rqu_et_sq_ave, rqu_et_sq_ave_sig_mean,rqu_et_sq_ave_intra,rqu_et_sq_ave_sig_mean_intra,rqu_et_sq_ave_inter,rqu_et_sq_ave_sig_mean_inter  ) )

        f_dG_r.write(" %f %d %16.12f %16.12f %16.12f  %16.12f  %16.12f  %16.12f \n"%(r_ij,N_cnts,dG_ave,dG_ave_sig_mean,dG_ave_intra,dG_ave_sig_mean_intra,dG_ave_inter,dG_ave_sig_mean_inter ))
        f_I_r.write(" %f %d %16.12f %16.12f %16.12f  %16.12f  %16.12f  %16.12f  \n"%(r_ij,N_cnts,I_ave,I_ave_sig_mean,I_ave_intra,I_ave_sig_mean_intra,I_ave_inter,I_ave_sig_mean_inter))
        f_lambda_r.write(" %f %d %16.12f %16.12f %16.12f  %16.12f  %16.12f  %16.12f  \n"%(r_ij,N_cnts,lambda_ave,lambda_ave_sig_mean,lambda_ave_intra,lambda_ave_sig_mean_intra,lambda_ave_inter,lambda_ave_sig_mean_inter))

        
        #f_VAB.write( "\n %f %d %16.12f  %16.12f  " % (r_ij , r_bin_cnt[b_indx], et_average, et_sig_mean ) )
        #  f_gv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , r_bin_cnt[b_indx], rhog_v, rhog_v_sig_mean, rhog_v_intra, rhog_v_sig_mean_intra, rhog_v_sig_mean_inter, rhog_v_sig_mean_inter  ) )

        # f_gv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , r_bin_cnt[b_indx], rhog_v, rhog_v_sig_mean, rhog_v_intra, rhog_v_sig_mean_intra, rhog_v_sig_mean_inter, rhog_v_sig_mean_inter  ) )


        # f_rqugv.write( "\n %f %d %16.12f  %16.12f " % (r_ij , r_bin_cnt[b_indx], rqu_rhog_v, rqu_rhog_v_sig_mean ) )

        #f_rqugv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , r_bin_cnt[b_indx], rqu_rhog_v, rqu_rhog_v_sig_mean, fourpi_rsqr_rhog_v_intra, fourpi_rsqr_rhog_v_sig_mean_intra, fourpi_rsqr_rhog_v_inter, fourpi_rsqr_rhog_v_sig_mean_inter  ) )

        
        # f_rsqgv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , r_bin_cnt[b_indx], fourpi_rsqr_rhog_v, fourpi_rsqr_rhog_v_sig_mean, fourpi_rsqr_rhog_v_intra, fourpi_rsqr_rhog_v_sig_mean_intra, fourpi_rsqr_rhog_v_inter, fourpi_rsqr_rhog_v_sig_mean_inter  ) )
        # f_dVgv.write( "\n %f %d %16.12f  %16.12f %16.12f  %16.12f  %16.12f  %16.12f " % (r_ij , r_bin_cnt[b_indx], dVg_v, dVg_v_sig_mean, fourpi_rsqr_rhog_v_intra, fourpi_rsqr_rhog_v_sig_mean_intra, fourpi_rsqr_rhog_v_inter, fourpi_rsqr_rhog_v_sig_mean_inter  ) )


        for dG_indx in range(dG_n_bins) :

		dG_ij = dG_indx*dG_bin_size + dG_floor # + r_bin_size*0.5
                # print b_indx,dG_indx ,dG_ij, dG_r_dG_bins[b_indx,dG_indx],N_cnts
		prob_dG = dG_r_dG_bins[b_indx,dG_indx]/float(N_cnts*2)
		f_dG_r3D.write( " %f %f %f \n"%(r_ij,dG_ij,prob_dG))
		
	f_dG_r3D.write( "\n" )

		 
	if( debug ):
	    if( b_indx == 50 ):
		print "N_cnts",r_bin_cnt[b_indx]
		print " et_sum ",et_sum
		print " et_average ",et_average
		print " et_sq ",et_sq
		print " et_sq_ave ",et_sq_ave
		print " et_qu_ave ",et_qu_ave
		print " et_sq_sig_mean ",et_sq_sig_mean
		print " r_out, r_in ",r_out, r_in
		print " dr_vol ",dr_vol
		print " group_cnt ",group_cnt
		print " dVg_v ",dVg_v
		print " dVg_v_sig_mean ",dVg_v_sig_mean
		print " rhog_v ",rhog_v
		print " rhog_v_sig_mean ",rhog_v_sig_mean
		print " fourpi_rsqr_rhog_v ",fourpi_rsqr_rhog_v
		print " fourpi_rsqr_rhog_v_sig_mean ",fourpi_rsqr_rhog_v_sig_mean
		
		print et_sq_ave, et_sq_sigma
		print rhog_v, rhog_v_sig_mean
		print fourpi_rsqr_rhog_v, fourpi_rsqr_rhog_v_sig_mean
		print dVg_v, dVg_v_sig_mean
		sys.exit(" debug mode! ")

    f_gr.close()
    f_VABsq.close()
    f_rquVABsq.close()
    f_dG_r.close()
    f_I_r.close()
    f_lambda_r.close()
    f_dG_r3D.close()
 
    print " Average I_rp_list %f eV "%(np.average(I_rp_list))
    print " Average I_pr_list %f eV "%(np.average(I_pr_list))

    l_sq = SUM_rsq_vsq / SUM_vsq
    print " r^2 V^2 / V^2 ",l_sq," l ",np.sqrt(l_sq)
    l_sq = SUM_rqu_vsq / SUM_rsq_vsq
    print " r^4 V^2 / r^ 2 V^2 ",l_sq," l ",np.sqrt(l_sq)
    
    N_pairs = sum(r_bin_cnt )
    N_Aq = N_pairs/volume_i
    
    print "   total pairs ",N_pairs,N_Aq,"  x10^-{21} "
    print "   Total: "
    print "     sum rho g_v ",SUM_rhog_v
    print "     sum r^2 rho g_v ",SUM_rsq_rhog_v
    print "     sum r^4 rho g_v ",SUM_rqd_rhog_v
    ave_l_sq = SUM_rqd_rhog_v/SUM_rsq_rhog_v
    ave_l_sq_sig_mean = np.sqrt( SUM_rqd_rhog_v_sig_mean) / np.sqrt( SUM_rsq_rhog_v_sig_mean )
    
    ave_l = np.sqrt(ave_l_sq )
    ave_l_sig_mean = np.sqrt(ave_l_sq_sig_mean )


    print "        shit snacks   <l^2> ",
    
    print "     <l^2> ",ave_l_sq, "+-",ave_l_sq_sig_mean,"     l ",ave_l, "+-",ave_l_sig_mean
    
    #
    # Intra mol 
    #
    
    N_pairs_intra = sum(r_bin_cnt_intra )
    N_Aq_intra = N_pairs_intra/volume_i*100
    
    print "   Intra molecular pair : ", N_pairs_intra,N_Aq_intra,"  x10^-{19}/cm^3 ",100.0*N_pairs_intra/N_pairs
    print "     sum rho g_v ",SUM_rhog_v_intra
    print "     sum r^2 rho g_v ",SUM_rsq_rhog_v_intra
    print "     sum r^4 rho g_v ",SUM_rqd_rhog_v_intra
    ave_l_sq = SUM_rqd_rhog_v_intra/SUM_rsq_rhog_v_intra
    
    ave_l_sq_sig_mean = np.sqrt( SUM_rqd_rhog_v_intra_sig_mean) / np.sqrt( SUM_rsq_rhog_v_intra_sig_mean )
    
    ave_l = np.sqrt(ave_l_sq )
    ave_l_sig_mean = np.sqrt(ave_l_sq_sig_mean )
    print "    Intra-molecular <l^2> ",ave_l_sq, "+-",ave_l_sq_sig_mean,"     l ",ave_l, "+-",ave_l_sig_mean
        
    #
    # Inter mol 
    #
    
    N_pairs_inter = sum(r_bin_cnt_inter )
    N_Aq_inter = N_pairs_inter/volume_i*100
    
    print "   Inter molecular pair : ", N_pairs_inter,N_Aq_inter,"  x10^-{19}/cm^3  ",100.0*N_pairs_inter/N_pairs
    print "     sum rho g_v ",SUM_rhog_v_inter
    print "     sum r^2 rho g_v ",SUM_rsq_rhog_v_inter
    print "     sum r^4 rho g_v ",SUM_rqd_rhog_v_inter
    ave_l_sq = SUM_rqd_rhog_v_inter/SUM_rsq_rhog_v_inter
    ave_l_sq_sig_mean = np.sqrt( SUM_rqd_rhog_v_inter_sig_mean) / np.sqrt( SUM_rsq_rhog_v_inter_sig_mean )
    
    ave_l = np.sqrt(ave_l_sq )
    ave_l_sig_mean = np.sqrt(ave_l_sq_sig_mean )
    print "    Inter-molecular <l^2> ",ave_l_sq, "+-",ave_l_sq_sig_mean,"     l ",ave_l, "+-",ave_l_sig_mean
        
    
if __name__=="__main__":
    main()
