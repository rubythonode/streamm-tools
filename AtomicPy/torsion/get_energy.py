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
    
    (options, args) = parser.parse_args()

    return options, args



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
