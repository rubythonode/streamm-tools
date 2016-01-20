#! /usr/bin/env python
# Make input files for torsional potential energy surface 

# Dr. Travis Kemper
# NREL
# 12/09/2013
# travis.kemper@nrel.gov

HtoeV = 27.211385

def get_options():
    import os, os.path
    from optparse import OptionParser
    usage = "usage: %prog [options] [input_files] \n"
    usage = usage + "Input files \n"
    usage = usage + "  specify the destination name of an option followed by the value"
    parser = OptionParser(usage=usage)
    

    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    
    parser.add_option("-r","--recalc", dest="recalc",action="store_true", default=False,help=" Rerun calculation even if finished calculation has been found ")

    # json files to act on
    parser.add_option("-j","--json", dest="json", default="",type="string",help=" json files to act on")

    (options, args) = parser.parse_args()
    
    return options, args


def pars_zmatrix( calc_id, job_name ):
    from string import replace
    # atomicpy
    import gaussian, top
    
    log_name = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.log" )
    com_name = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.com" )
    # Get lines of log file     
    f = open(log_name,'r')
    Lines = f.readlines()
    f.close()

    #Parse fchk
    fchk_file = "%s/%s%s" % ( calc_id ,job_name, "-ZMATOPT.fchk" )

    print " fchk_file" , fchk_file
    NA, ELN, R, TOTAL_ENERGY, Q_ESP  = gaussian.parse_fchk( fchk_file )
    
    # Parse log file 
    NBLIST,NBINDEX = top.build_covnablist(ELN,R)
    BONDS  = top.nblist_bonds(NA,NBLIST, NBINDEX)
    
    # Find rings 
    RINGLIST, RINGINDEX , RING_NUMB = top.find_rings(ELN,NBLIST,NBINDEX)
    RING_CONNECT  = top.find_conections(ELN,NBLIST,NBINDEX,RINGINDEX , RING_NUMB) 

    zmatrix = gaussian.com_zmatrix(com_name)    
    
    DIH_ID, DIH_VAL, DIH_ATOMS = gaussian.get_dih_id( zmatrix)
    
    return ( RING_CONNECT, RING_NUMB, DIH_ID, DIH_VAL, DIH_ATOMS, zmatrix )


def main():
    import string, os , sys 
    # atomicpy
    import jsonapy
    import gaussian, elements, xmol , file_io , cluster , top 
    from string import replace
    
    #
    # Set some defaults 
    #
    default_method = 'b3lyp'
    default_basis = '6-31G**'
    
    options, args = get_options()
    
    work_dir = os.getcwd()

	
    json_files = options.json.split(',')
    print json_files
    if( len(json_files) > 0 ):
	# Read index files from args
	for json_file in json_files:
	    # Get lines of index file
		
	    
	    # Verbose output
	    if( options.verbose ):
		print "The molecules specified in json file ",options.json," will be read in "
    
	    json_data,json_success = jsonapy.read_jsondata(json_file)
	    if(  json_success ):
		#
		mol_dir,tag,n_units,accuracy,method,basis,acceptors,acceptor_substituents,donors,donor_substituents,terminals,terminal_substituents,spacers,spacer_substituents,metadata_found = jsonapy.read_meta(json_data)
		#
		# Need meta data to proceed 
		#      		    
		if( metadata_found ):
		    json_atomicdata = 0
		    fchk_atomicdata = 0
		    
		    if( options.verbose ):
			print " Meta data found will use specified method and basis unless others are specified in the options "
		    #
		    # Construct file names 
		    #
		    #short_name = "acc%d_%s_n%d" % (accuracy, tag, number )
		    job_name = "acc%d_%s_n%d" % (accuracy, tag, n_units )
		    struct_dir = "%s/%s/" % (mol_dir, tag )
		    calc_id = "%s/%s%s" % (struct_dir, job_name , "-ZMATOPT" )
		    #
		    # 
		    #
		    zmat_fchk = "%s/%s%s" % ( calc_id ,job_name,"-ZMATOPT.fchk" )
		    print " Checking for complete zmatrix optimiztion ",zmat_fchk
                    zmat_finished = file_io.file_exists( zmat_fchk )
		    
		    if(  zmat_finished ):
                        
                        if( options.verbose ):
                            print " Reading atomic data from ",zmat_fchk
                        NA, ELN, R, TOTAL_ENERGY , Q_ESP , N_ELECTRONS,NMOS,EIGN = gaussian.parse_fchk_HL( zmat_fchk )
                        en_line = "%d  %f %f \n "%( n_units,EIGN[N_ELECTRONS-1],EIGN[N_ELECTRONS])

			hl_qm = struct_dir +'/' +tag + '_HL.dat'
                        if( options.verbose ):
                            print " Writing data file ",hl_qm
                            
                        if( options.recalc ):
                            hl_out = open(hl_qm,'w')
                            hl_out.write( '# units , HOMO (eV) , LUMO (eV) \n')
                            hl_out.close()

                        hl_out = open(hl_qm,'a')
                        hl_out.write(en_line)
                        hl_out.close()
                        

                            
                
	
	

		    else:
			print "  Zmatrix optimization has not finished for ",job_name
			    						    
		
		    
    else:
	print " No json files specified "
        
    
if __name__=="__main__":
    main() 

