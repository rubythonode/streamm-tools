


const_avo = 6.02214129 # x10^23 mol^-1 http://physics.nist.gov/cgi-bin/cuu/Value?na

def get_options():
    import os, os.path
    from optparse import OptionParser
 
    usage = "usage: %prog [options] \n"
    parser = OptionParser(usage=usage)
    
    parser.add_option("-v","--verbose", dest="verbose", default=False,action="store_true", help="Verbose output ")
    parser.add_option("--in_data", dest="in_data", type="string", default="", help="Input lammps structure file (.data) ")
    parser.add_option("--in_lammpsxyz", dest="in_lammpsxyz", type="string", default="", help="Input lammps xyz file with atoms listed as atom type numbers")

    parser.add_option("--out_xyz", dest="out_xyz", type="string", default="", help=" Output xyz file ")
    parser.add_option("--out_data", dest="out_data",type="string",default="",help=" Output Lammps data file ")
    #
    # xmol output 
    #
    parser.add_option("--out_xmol", dest="out_xmol", type="string", default="", help=" Output xmol file ")
    parser.add_option("--frame_o", dest="frame_o", type=int, default=0, help=" Initial frame to read")
    parser.add_option("--frame_f", dest="frame_f", type=int, default=0, help=" Final frame to read")
    parser.add_option("--frame_step", dest="frame_step", type=int, default=1, help=" Read every nth frame ")
    
    (options, args) = parser.parse_args()
        
    return options, args
   
def main():
    #
    # Caclulate rdf
    #
    
    import os, sys, numpy , math , random, json
    import datetime
    import time    
    import gromacs, elements, xmol, prop, file_io, groups,lammps , top, jsonapy

    debug = 0
    
    # Load information onto all processors 
    #
    options, args = get_options()
    prop_dim = 3



    lammpsxyz_F = open(options.in_lammpsxyz , 'r' )
    lammpsxyz_lines = lammpsxyz_F.readlines()
    lammpsxyz_F.close()        


    NA_sys=int(lammpsxyz_lines[0] )
    options.frame_f  = int( float( len(lammpsxyz_lines) )/ float( NA_sys + 2) )

    print " Analyzing ",options.in_lammpsxyz," with ",NA_sys," atoms and ",options.frame_f ," frames "

    frame_cnt = 0
    for frame_i in range(options.frame_o,options.frame_f,options.frame_step):
        line_cnt = frame_i*(NA_sys + 2 ) 
        frame_cnt += 1 
        if( len(options.out_xmol) ):
            str_file.write( str(NA_sys) + "\n" )
            
        line_cnt += 1
                
        if( options.verbose ):
            print " reading frame ",frame_i," starting at line ",line_cnt-1," with comment ",lammpsxyz_lines[line_cnt] 
                
        if( len(options.out_xmol) ):
            str_file.write( lammpsxyz_lines[line_cnt] )

        TYPES_f = []
        for atom_i in range(NA_sys ):
            line_cnt += 1
                    
            if( line_cnt > len(lammpsxyz_lines)-1):
                print " frame is missing some atoms ",atom_i," not found "
                # sys.exit("read in e)
                        
            col =  lammpsxyz_lines[line_cnt].split()
            if( len(col) >= 4 ):

                type_i = int(col[0]) 
                TYPES_f.append(type_i)

        if( frame_cnt > 1 ):
            print " check atom types against previous frame "
            for atom_i in range(len(TYPES_f) ):
                if(  TYPES_f[atom_i] !=  TYPES_pf[atom_i] ):
                    print " atom ",atom_i," has changed from ", TYPES_f[atom_i]," to ",  TYPES_pf[atom_i]
                    sys.exit(" bad file ")
                    
            
        TYPES_pf = []
        for atom_i in range(len(TYPES_f) ):
            TYPES_pf.append( TYPES_f[atom_i] )


if __name__=="__main__":
    main()
   
