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

    parser.add_option("-i","--input", dest="input", type="string", default="et.dat", help="Input data file ")

    # Gromacs
    parser.add_option("-o","--output", dest="output", type="string", default="out.hist", help="Output histogram file ")
    
    (options, args) = parser.parse_args()

    return options, args

def main():
    import os, sys, numpy , math 
    import file_io
    #import statistics
    
    options, args = get_options()
    
    print " Read in file " , options.input
    
    data_file =  options.input
    if( file_io.file_exists(data_file)):
        dat = open(data_file,"r")
        Lines = dat.readlines()
        dat.close()

    bin_size = 1.0 #(angstroms)
    v_min = 0.0
    v_max = 100.0
    v_floor = int(v_min/bin_size)*bin_size
    val_range = v_max - v_floor
    n_bins = int(val_range/bin_size) + 1

    bins = numpy.zeros(n_bins)
    val_sq_bins  = numpy.zeros(n_bins)
    bin_cnt = numpy.zeros(n_bins)

    
    col_numb = 1
    
    for line in Lines:
        col = line.split()
        if( col[0][0] != "#" and col[0][0] != "@"  ):
            t_i = float( col[0] )
            val_i = float( col[col_numb] )
	    
            bin_v = int( val_i/ bin_size )
            if( bin_v < n_bins ):
                bins[bin_v] += val_i
                val_sq = val_i*val_i 
                val_sq_bins[bin_v] += val_sq
                bin_cnt[bin_v] += 1
		
		
    val_sum = sum(  bins )
    val_sq_sum = sum(  val_sq_bins )

    cnt_sum = sum(bin_cnt)
    ave = val_sum/float(cnt_sum) 
    
    print cnt_sum, val_sum , val_sq_sum
    
    val_stdv = numpy.sqrt( val_sq_sum/cnt_sum - ( val_sum/cnt_sum)**2 ) 
    
    f_data = open(options.output,"w")
    f_data.write("# Average %f " % ave )
    f_data.write("# standard deviation  %f " %  val_stdv )
    for b_indx in range(n_bins):
        y_i = b_indx*bin_size + v_floor + bin_size*0.5 
	
        if(  bin_cnt[b_indx] > 0 ):
            b_sum =  bins[b_indx] 
            b_prop = b_sum/float( val_sum )
	    b_average = b_sum/float( bin_cnt[b_indx] )

            val_sq = val_sq_bins[b_indx]
            val_sq_ave = numpy.sqrt(  val_sq_sum/float( bin_cnt[b_indx] ) - (b_average)**2 ) 
            
	    
                
            #ket_value = ket(b_average)q
        else:
            b_sum = 0.0 
            b_prop = 0.0 
            val_sq  = 0.0 
            val_sq_ave = 0.0
	    
	f_data.write( "\n %f %f %f " % (y_i, b_prop , val_sq_ave ) )
     
    f_data.close()
    
    print " Average ",val_sum/float(cnt_sum),"  standard deviation  ",val_stdv

if __name__=="__main__":
    main()
