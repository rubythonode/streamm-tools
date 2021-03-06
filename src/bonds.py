"""
Class data structures for 2, 3, 4 point groupings of Particle objects
"""
import copy, sys

class Bond:
    """
    Data structure for describing any 2-point associatiaon of Particle-s
    """

    def __init__(self, pgid1=0, pgid2=0, length=0.0, type="blank"):
        """
        Constructor for a general bond. Checks for types in arguments
        and throws a TypeError when appropriate

        Args:
            pgid1   (int)   GlobalID of Particle object in bond
            pgid2   (int)   GlobalID of Particle object in bond
            length (float) Cartesian length of bond
            type   (str)   Charge value in units of [e]
        """
        
        if isinstance(pgid1, int):
            self.pgid1 = pgid1
        else:
            print "1st arg should be int"
            raise TypeError

        if isinstance(pgid2, int):
            self.pgid2 = pgid2
        else:
            print "2nd arg should be int type"
            raise TypeError

        if isinstance(length, float):
            self.length = length
        else:
            print "3rd arg should be float value"
            raise TypeError

        if isinstance(type, str):
            self.type = type
        else:
            print "4th arg should be string value"
            raise TypeError

        self.lmpindx = 0 
        self.g_indx = 0 


    def __del__(self):
        """
        Destructor, clears object memory
        """
        del self.pgid1
        del self.pgid2
        del self.length
        del self.type
        del self.lmpindx
        del self.g_indx 

    def __contains__(self, pgid):
        """
        'Magic' method implementing 'in' keyword.

        Args:
            pgid (int) Particle GID to check against 'held' IDs
        """
        if ( (pgid == self.pgid1) or (pgid == self.pgid2) ):
            return True
        else:
            return False

    def set_lmpindx(self,lmpindx):
        """
        Set bond type index for lammps
        """
        self.lmpindx = lmpindx
        
        
    def get_lmpindx(self):
        """
        Return bond type index for lammps
        """
        return self.lmpindx


    def set_g_indx(self,g_indx):
        """
        Set bond type index for gromacs 
        """
        self.g_indx = g_indx
        
        
    def get_g_indx(self):
        """
        Return bond type index for gromacs
        """
        return self.g_indx

    def __str__(self):
        """
        'Magic' method for printng contents
        """
        return " %s - %s    %s "%(self.pgid1,self.pgid2,self.type )

class BondContainer:
    """
    Main data structure for holding Bond objects. Map of global
    bond ID (integer) to Bond object instances
    """

    def __init__(self, idList=[], verbose=False):
        """
        Constructor: sets up a dictionary for indexing 'Bond' objects

        Args:
            idList (list): of bond IDs. If empty then ID starts at 1.
                If not empty then ID's (keys) are inititalized with Bond objects
            verbose (bool): flag for printing status/debug info        
        """
        self.verbose=verbose
        self.bonds=dict()                          # Creates empty dict struc
        self.bonds={key: Bond() for key in idList} # Creates empty Bond objs
                                                   #   if idList not empty

        if len(idList) == 0:         # If list not set in constructor arg
            self.maxgid=0            # default=0 if idList empty
        else:                        #
            self.maxgid=max(idList)  # take max in list for maxgid


    def __del__(self):
        """
        Destructor, clears dictionary memory
        """
        if self.verbose:
            print "Cleaning particle container"
        del self.bonds
        del self.maxgid


    def clear(self):
        """
        Clears bonds out of BondContainer
        """
        self.maxgid = 0
        self.bonds=dict()                          # Creates empty dict struc


    def __len__(self):
        """
        'Magic' method for returning size of container
        """
        return len(self.bonds)

    def __str__(self):
        """
        'Magic' method for printng contents
        """

        bondStr="\n Contains bond objects: \n"
        for gid in self.bonds:
            bondStr += " %d : %s \n"%(gid,self.bonds[gid] )

        return bondStr


    def keys(self):
        """
        Return list of all ptcl IDs (keys) currently in container
        """
        keyList = self.bonds.keys()
        return keyList


    def __setitem__(self, gid, bond):
        """
        'Magic' method implementing obj[]=value operator
        Performs deep copy of value so container is managing memory
        """
        if gid in self.bonds.keys():
            self.bonds[gid]=copy.deepcopy(bond)
        else:
            print "Cannot add bond object to non-existent ID"
            sys.exit(3) 


    def __getitem__(self, gid):
        """
        'Magic' method implementing obj[] operator
        Operations on returned elements change container
        """
        return self.bonds[gid]


    def __delitem__(self, gid):
        """
        'Magic' method implementing del obj[] operator
        """
        del self.bonds[gid]


    def __iter__(self):
        """
        'Magic' method implementing (for x in 'this')....
        """
        # return iter(self.bonds)
        return self.bonds.iteritems()


    def __call__(self, idSubList=None):
        """
        Callable magic method. Returns iterator to subset bonds dictionary

        Args:
             idSubList (list) list of pid-s of particle objects to be returned
             
        Returns: iterator to subset of particle dictionary
        """

        subGroupDct = dict()
        
        if idSubList != None:
            for gid, bondObj in self.bonds.iteritems():
                if gid in idSubList:
                    subGroupDct[gid] = bondObj
            return subGroupDct.iteritems()

        else:
            print "Callable BondContainer requires a list of subgroup bond IDs"
            sys.exit(3)


    def __contains__(self, gid):
        """
        'Magic' method implementing in keyword (key in obj')....
        """
        return gid in self.bonds


    def hasBond(self, bondList):
        """
        Check the ptcl IDs in bondList for any bond in container that is similar
        eg bond 1-2 is same as bond 2-1

        Args: (list) ptcl IDs defining bond to search for
        
        Returns: (bool) is bond in container
        """

        for gid, bondObj in self.bonds.iteritems():
            p1 = bondObj.pgid1
            p2 = bondObj.pgid2
            if ((p1 in bondList) and (p2 in bondList)):
                return True
            
        return False


    def deletepid(self, pid):
        """
        Delete bond with a partile id 
        """
        del_id = -1 
        for gid, bondObj in self.bonds.iteritems():
            p1 = bondObj.pgid1
            p2 = bondObj.pgid2
            if ( p1 == pid or p2 == pid ):
                #del bondObj
                del_id = gid
                # print " deleting ",gid,p1,p2
        if( del_id > 0 ):
            self.bonds.pop(del_id)
        else:
            error_line = "!!! Pid {} not found in bonds !!!".format(pid)
            print error_line
        return
    

    def __iadd__(self, other):
        """
        'Magic' method to implement the '+=' operator
        
        Compare global IDs of bonds and reassign globalIDs for bond
        container using the max ID between the two lists

        Note: for now this reassigns ID always
        """

        keys1 = self.bonds.keys()         # global IDs in this object
        keys2 = other.bonds.keys()        # global IDs in object being added
        bothkeys = keys1 + keys2          # List of all keys

        if len(bothkeys) > 0:                 # If keys not empty... proceed
            self.maxgid = max(keys1 + keys2)  # find max globalID in keys, set this object maxID

        for ptclkey2 in other.bonds:
            self.put(other.bonds[ptclkey2])

        return self


    def put(self, bond):
        """
        Append 'Bond' object to this container. Updates globalID for container
        by incrementing the maxgid member

        Args:
            ptcl (Particle) correctly initialized Particle object

        NOTE:
            (1) One can imagine extra conditions on bonds inserted
            (2) This could check for uniqueness of all globalID's and throw error for copies
        """
        
        if isinstance(bond, Bond):
            self.maxgid += 1
            self.bonds[self.maxgid] = copy.deepcopy(bond)
        else:
            print "Attempting to add non-Bond type to container"
            raise TypeError


    def replacePtclIDs(self, idFromTo):
        """
        Replace ptclIDs given a dictionary of ID changes # eg {1:3, 3:5, 2:20...}
                
        Args:
            idFromTo (dict) map of ID changes
        """

        fromIDs = idFromTo.keys()
                
        for gid in self.bonds:
            
            bond = self.bonds[gid]   # Bond object
            pgid1 = bond.pgid1       # ptcl1 in bond
            pgid2 = bond.pgid2       # ptcl2 in bond

            if pgid1 in fromIDs:
                toID = idFromTo[pgid1]
                bond.pgid1 = toID

            if pgid2 in fromIDs:
                toID = idFromTo[pgid2]
                bond.pgid2 = toID




    def getTypeInfoDict(self):
        """
        Return a map of type to typeIndex
        Method assigns a type index and checkes for consistency

        Returns:
            dictionary
        """

        # Look for types and get unique list
        typeList = list()
        for gid, bondObj in self.bonds.iteritems():
            bondType = bondObj.type
            typeList.append(bondType)
        typeList = list(set(typeList))
        
        # Generate list of unique type for keys to initialize dictionary
        typeIndexDict = {key:index+1 for index, key in enumerate(typeList)} # dict for typeIndex
        if self.verbose:
            print "Unique types = ", typeList
            print "typeIndexDict = ", typeIndexDict

        # Pack yypeIndex to new dictionary
        typeInfoDict = dict()
        for key in typeIndexDict:
            index  = typeIndexDict[key]
            typeInfoDict[key] = index

        return typeInfoDict
