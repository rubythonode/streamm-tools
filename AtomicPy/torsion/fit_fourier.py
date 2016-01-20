#! /usr/bin/env python

# run ab initio torsional potential 

# Dr. Travis Kemper
# NREL
# 12/09/2013
# travis.kemper@nrel.gov

EVTOKCAL = 23.0605

def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    usage = usage + "Input files \n"
    usage = usage + "  specify the destination name of an option followed by the value"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    # Cluster options
    parser.add_option("--host", dest="host",type="string",default="macbook",help=" name of machine  ")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")

    parser.add_option("--gnuplot", dest="gnuplot", default=False, action="store_true", help="Output gnuplot file ")
    
    parser.add_option("--tor_paramo", type="float", dest="tor_paramo", default=0.09, help="Value to initialize torsional potential coefficients ")
    parser.add_option("--out_ftor", dest="out_ftor", type="string", default="tor_fit.dat", help="Output fitted torsional potential data ")
    parser.add_option("--plot_tor", dest="plot_tor", type="string", default="tor_fit.pdf", help="Output fitted torsional potential graph ")
    parser.add_option("--plot_tor_comp", dest="plot_tor_comp", type="string", default="tor_comp.pdf", help="Output for components of torsional potential graph ")

    parser.add_option("--qm_tor", dest="qm_tor", type="string", default="",  help=" Data file of quantum target values ")
    parser.add_option("--ff_tor", dest="ff_tor", type="string", default="",  help=" Data file of force field values ")
    
    (options, args) = parser.parse_args()

    return options, args


def plot_fourierfit(Fourier_coef, angle_list, targets,fourier_fit_plt, tag, n_units, d_pos ,rmse_four ):
    import matplotlib.pyplot as plt
    import math, sys 
    
    debug = 0
    
    target_kcal = []
    fit_en_kcal = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        target_kcal.append( t_en*EVTOKCAL)
        ang_val = angle_list[cent_indx] #[0]
        theta = math.radians( ang_val )
        tor_en = FourierSeries(Fourier_coef,theta)
        fit_en_kcal.append( tor_en*EVTOKCAL )
    
    # plt.ylabel('Energy (kcal/mol)')
    tittle_lab = tag + " n=" + str(n_units) 
    
    plt.suptitle(tittle_lab)
    plt.ylabel('Energy (kcal/mol)')
    plt.xlabel('dihedral angle (deg)')
    
    rmse_fourkcal = rmse_four*EVTOKCAL
    rmse_fourkcal_r = round(rmse_fourkcal , 4 )
    
    t_lab = "target p="+str(d_pos) + "\n RMSE = " + str(rmse_fourkcal_r) + " kcal/mol "
    f_lab = "fit p="+str(d_pos)
    plt.plot( angle_list,target_kcal,'k-', label=t_lab )
    plt.plot( angle_list,fit_en_kcal,"rx", label=f_lab)
    plt.legend(loc=(0.67,0.72),prop={'size':10})
    
    #plt.show()
    
    plt.savefig(fourier_fit_plt,format='pdf')
    plt.close()


def qm_analysis(  dih_qm ):
    import file_io
    
    success = 0 
    
    qm_min =  1e16 
    qm_max = -1e16
    ang_min = "nan"
    ang_max = "nan"
    
    #ang_min 
    
    angle_list = []
    energy_list = []
    
    
    if( file_io.file_exists(dih_qm) ):
	    
	f_qm = open(dih_qm ,"r")
	f_qm_Lines = f_qm.readlines()
	f_qm.close()
	
	
	for f_qm_line in f_qm_Lines:
	    f_qm_col = f_qm_line.split()
	
	    if( len(f_qm_col) >= 3 and f_qm_col[0] != "#" ):
		qm_angle = float( f_qm_col[1] )
		qm_en = float( f_qm_col[2] )
                
                angle_list.append(qm_angle)
                energy_list.append(qm_en)
                
		success = 1 
	    
    return (success, angle_list,energy_list )

def FourierSeries(p,theta_rad):
    import math
    # theta_rad dihedral angle in radians
    
    debug = 0
    
    fourier_sum = 0.0
    
    for n in range(len(p)):
        fourier_sum +=  p[n]*( math.cos( float(n)*theta_rad  ) )
        if( debug ):
            print n,p[n],theta_rad, math.cos( float(n)*theta_rad ), p[n]*( math.cos( float(n)*theta_rad  ) )
    
    return fourier_sum

def FourierSeries_comp(p,theta_rad):
    import math
    # theta_rad dihedral angle in radians
    
    debug = 0
    
    fourier_comp = []
    
    for n in range(len(p)):
        fourier_comp.append(  p[n]*( math.cos( float(n)*theta_rad  ) ) )
    
    return fourier_comp

def residuals_FS(Fourier_coef, angle_list, targets,wt_angle,wt_coef ):
    import math, sys 
    
    debug = 0
    
    # Round parameters to ~0.01 kca/mol
    for param_i in range( len(Fourier_coef) ):
        rounded_param = round(Fourier_coef[param_i] ,8 )
        Fourier_coef[param_i]   = rounded_param 
    
    resid = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        ang_val = angle_list[cent_indx] #[0]
        wt = wt_angle[cent_indx]

        theta = math.radians( ang_val )
        tor_en = FourierSeries(Fourier_coef,theta)
        
        delta_en =  t_en - tor_en
        sq_delta = wt*( delta_en*delta_en)
        
        resid.append(sq_delta)
        
    return resid


def rmse_FS(Fourier_coef, angle_list, targets ):
    import math, sys 
    
    debug = 0
    
    resid = []
    fit_en = []
    for cent_indx in range(len(targets)):
        t_en = targets[cent_indx]
        ang_val = angle_list[cent_indx] #[0]
        theta = math.radians( ang_val )
        tor_en = FourierSeries(Fourier_coef,theta)
        fit_en.append( tor_en )
        delta_en =  t_en - tor_en
        sq_delta =  delta_en*delta_en
        
        resid.append(sq_delta)
        
    rmsa = sum( resid)/float( len(targets) )
    rmse = math.sqrt( rmsa )
    
    return rmse

def set_weights(angle_list,eneryg_list,min_indx,max_indx):
    import sys
    import numpy as np
    
    debug = 0 
        
    wt_o = 10.0 
    wt_max = 1000.0 
    wt_min = 1000.0 
    wt_min_pm = 30.0
    wt_min_pm2 = 20.0 

    wt_coef = 0.0
    wt_angle = []

    n_anlges =  len(angle_list)
    
    wt_angle = [] #np.ones(n_anlges+1)
    for cent_indx in range(n_anlges+1):
        wt_angle.append(wt_o)
    
    #h = angle_list[1] - angle_list[0]
    #success,  tor_en,min_indx,max_indx,trans_list, k_list = prop.calc_d2(  qm_en_s , h)


    for indx in range( len(min_indx) ):
        
        cent_indx =  min_indx[indx]
        wt_angle[cent_indx] = wt_min
        if( cent_indx < len(angle_list)  ):  wt_angle[cent_indx + 1 ] = wt_min_pm
        if( cent_indx < len(angle_list) - 1 ):  wt_angle[cent_indx + 2 ] = wt_min_pm2
        if( cent_indx > 0 ):  wt_angle[cent_indx - 1 ] = wt_min_pm
        if( cent_indx > 1 ):  wt_angle[cent_indx - 2 ] = wt_min_pm2
        
        if(debug): print "  Setting min at ",angle_list[cent_indx]," to ",wt_min #," with dE = ",k_list[indx]
    
    for indx in range( len(max_indx) ):
        
        cent_indx =  max_indx[indx]
        wt_angle[cent_indx] = wt_max
        
    debug = 0
    if( debug ): 
        for cent_indx in range(len(angle_list)):
            print  angle_list[ cent_indx ], eneryg_list[ cent_indx ], wt_angle[cent_indx]
        sys.exit(" debug weights ")
        
    wt_coef = 1.0
    
    return ( wt_angle, wt_coef)

    

def main():
    import os, sys
    import jsonapy, prop
    import collections
    import math, numpy 
    from scipy import optimize
    
    debug = 0 
    n_max  = 14

    EVTOKCAL = 23.0605
    EVTOkJ = 96.4853
    
    options, args = get_options()

    qm_sufix = "_qm2"  # Need to read from qm_sufix_list[dih_indx]

    # Read index files from args
    for indx_file in args:
        # Get lines of index file   
        f = open(indx_file,'r')
        Lines = f.readlines()
        f.close()
	
	
	# Initialize lists
	style_lines = []
	plot_lines = []
	calc_i = 0 
	
        for line in Lines:
            col = line.split()
            if( len(col) >= 4 and col[0] != "#" ):
		
		
		struct_dir =  col[3]
		job_name = col[4]
		
		json_file =  struct_dir + "/" + job_name +".json"
		json_data,json_success = jsonapy.read_jsondata(json_file)
		if(  json_success ):
		    
		    mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)		
		    #
		    # Need meta data to proceed 
		    #      		    
		    if( metadata_found  ):
	    
			print "   Plotting ",struct_dir,job_name
			# Open output file
				
			qm_dih_id_list ,qm_cent_min_list ,qm_cent_max_list ,qm_cent_step_list,qm_a_k_list, qm_a_i_list, qm_a_j_list, qm_a_l_list,qmtor_found = jsonapy.read_qm_tor(json_data)
			if( qmtor_found ):
			    bar_cnt = [-1]*n_max*2 # numpy.zeros( [n_max*2])

			    for dih_indx in range( len(qm_dih_id_list) ):
				dih_id = qm_dih_id_list[dih_indx]
				cent_min = qm_cent_min_list[dih_indx]
				cent_max = qm_cent_max_list[dih_indx]
				cent_step = qm_cent_step_list[dih_indx]
				
                                
				dih_qm = struct_dir +'/' +job_name + '-' + dih_id + qm_sufix + '.dat'
				
				success, angle_list,energy_list = qm_analysis(  dih_qm )
		    

				bar_cnt[n_units] += 1 
                                    
				if( success ):
                                    print "  Data file ",dih_qm
                                    delta_param = True
                                    fit_iter = 0
                                    
                                    targets = []
                                    for angle_indx in range( len(energy_list) ):
                                        targets.append( energy_list[angle_indx] - min(energy_list) )
                        
                        
                                    fourier_order = 6
                                    Fourier_coef_fit = []
                                    param_list_pre = []
                                    for n in range(fourier_order+1):
                                        Fourier_coef_fit.append( 0.0 )
                                        param_list_pre.append(  0.0 )
        
                                    # Find second derivates for max/min
                                    h = angle_list[1] - angle_list[0]
                                    success,  tor_en,min_indx,max_indx,trans_list,trans_indxs, k_list = prop.calc_d2(  targets , h)
                                    
                                
                                    # Minimum 
                                    for indx in range( len(min_indx) ):
                                        cent_indx =  min_indx[indx]
                                        if( cent_indx > 0 ):
                                            print " Minimum target_en_s meV ",targets[cent_indx-1]*1000 , targets[cent_indx]*1000 ,  targets[cent_indx+1]*1000
        
                                    # Set wieghts
                                    wt_angle, wt_coef = set_weights(angle_list,targets,min_indx,max_indx)
                                    #
                                    
                                    print " fitting ",len(Fourier_coef_fit)," parameters to data ", len(targets)
                                    while delta_param:
                                        fit_iter +=  1
                                        
                                        Fourier_coef_fit,success = optimize.leastsq(residuals_FS,Fourier_coef_fit,args=(angle_list, targets,wt_angle,wt_coef ),epsfcn=0.0001)
                                    
                                        d_param = 0.0 
                                        for p_indx in range( len(Fourier_coef_fit)):
                                            d_pf = Fourier_coef_fit[p_indx] - param_list_pre[p_indx]
                                            d_param += numpy.sqrt(d_pf*d_pf )
                                        
                                        print "   fit_iter ",fit_iter," with delat_param ",d_param
                                        if( d_param <  0.0001 ):
                                            delta_param = False
                                        
                                        param_list_pre = []
                                        for p_indx in range( len(Fourier_coef_fit)):
                                            param_list_pre.append( Fourier_coef_fit[p_indx] )
                                                
                                    # if( options.verbose ):
    
                                    #
                                    # Print data file with fitted energies 
                                    #
                                    dih_fourier = struct_dir +'/' +job_name + '-' + dih_id  + '_fourier.dat'
                                    fourier_fit_plt = struct_dir +'/' +job_name + '-' + dih_id  + '_fourierfit.pdf'
                                    
                                    f_fourier  = open(dih_fourier,'w')
                                    f_fourier.write( "# ang_val,t_en,tor_en,sq_delta " )
                                    
                                    for cent_indx in range(len(targets)):
                                        t_en = targets[cent_indx]
                                        ang_val = angle_list[cent_indx] #[0]
                                        # wt = wt_angle[cent_indx]
                                        wt = 1.0 
                                
                                        theta = math.radians( ang_val )
                                        tor_en = FourierSeries(Fourier_coef_fit,theta)
                                        fourier_comp = FourierSeries_comp(Fourier_coef_fit,theta)
                                        
                                        delta_en =  t_en - tor_en
                                        sq_delta =  delta_en*delta_en
                                        
                                        #if( options.verbose ):
                                        #    print ang_val,t_en,tor_en,sq_delta
                                        
                                        f_fourier.write( "\n %f %f %f %f  %f %f %f %f  %f %f %f " %  (ang_val,t_en,tor_en,sq_delta,fourier_comp[0],fourier_comp[1],fourier_comp[2],fourier_comp[3],fourier_comp[4],fourier_comp[5],fourier_comp[6]))
                                        
                                    f_fourier.close()
                                    
                                    if( options.gnuplot ):
                                        print " plot \'",dih_fourier,"\' us 1:2 title \'target\'"
                                        print "  \'\' us 1:3 title \'fit\'"
                                    
                                    rmse_four = rmse_FS(Fourier_coef_fit, angle_list, targets )
                                    
                                    f_string = ""
                                    for n in range(fourier_order+1):
                                        f_string += str( Fourier_coef_fit[n] ) + " "
                                    f_string += str( rmse_four) + " "
                                    print " Fourier coef for ",dih_qm,f_string
                                    
				    d_o = 1*(n_units-2)
				    
				    dih_cnt =  int( bar_cnt[n_units] )
				    d_indx = n_max + d_o + dih_cnt*2
				    d_pos =   d_o - dih_cnt*2 
				    
                                    
                                    tag_fourier = struct_dir +'/' + job_name + '_fourier.dat'
                                    tag_file = open(tag_fourier,"a")
                                    tag_file.write(  "\n %d %d %s %s " %  (n_units,d_pos,dih_id,f_string))
                                    tag_file.close()
                                    
                                    if( d_pos == 0 or d_pos == 1 ):
                                            
                                        tag_fourier_cent = struct_dir +'/' + tag  + '_cent_fourier.dat'
                                        tag_file_cent = open(tag_fourier_cent,"a")
                                        tag_file_cent.write(  "\n %d %d %s %s " %  (n_units,d_pos,dih_id,f_string))
                                        tag_file_cent.close()
                                        
                                        print "rm ",tag_fourier," ",tag_fourier_cent
                                        
                                    #print " Plotting fit in ",fourier_fit_plt
                                    #plot_fourierfit(Fourier_coef_fit, angle_list, targets,fourier_fit_plt,tag, n_units,d_pos,rmse_four)
                                        
                                    # print "plot \'",tag_fourier_cent,"\' us 1:4,'' us 1:5, '' us 1:6,'' us 1:7,'' us 1:8, '' us 1:9, '' us 1:10"
                                    
if __name__=="__main__":
    main() 
