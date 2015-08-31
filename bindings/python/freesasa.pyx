"""FreeSASA Python interface.

This module provides a python interface to FreeSASA

Example
-------
A minimal program would be something like

    structure = freesasa.Structure("1abc.pdb")
    result = freesasa.calc(structure)
    print "SASA = %f" % result.totalArea()

See documentation of the classes and functions for how to customize behavior


The private methods _assign_ptr allow access to the C structs that
these classes are wrapping, which are necessary in the functions
calc() and classifyResults(). They are not intended to be used outside
of this module.
"""
from libc.stdio cimport FILE, fopen, fclose
from libc.stdlib cimport free, realloc
from libc.string cimport memcpy
from cpython cimport array
from cfreesasa cimport *

ShrakeRupley = 'ShrakeRupley'
"""string: Used to specify the algorithm by Shrake & Rupley"""

LeeRichards = 'LeeRichards'
"""string: Used to specify the algorithm by Lee & Richards"""

polar = FREESASA_POLAR

apolar = FREESASA_APOLAR

silent = FREESASA_V_SILENT
"""int: Suppress all warnings and errors (used by setVerbosity*())"""

normal = FREESASA_V_NORMAL
"""int: Normal verbosity (used by setVerbosity())"""

defaultParameters = {'algorithm' : ShrakeRupley, #freesasa_default_parameters.alg,
                     'probe-radius' : freesasa_default_parameters.probe_radius,
                     'n-points' :
                       freesasa_default_parameters.shrake_rupley_n_points,
                     'delta' :
                       freesasa_default_parameters.lee_richards_delta,
                     'n-threads' : freesasa_default_parameters.n_threads }
"""The default values for calculation parameters"""


cdef class Parameters:
      """
      Stores parameter values to be used by calculation.
      
      Wraps the C struct freesasa_parameters.
      """
      cdef freesasa_parameters _c_param
      def __init__(self,param=None):
            """
            Initializes Parameters object.
            
            Args:
               param (dict): optional argument to specify parameter-values, modeled after
                  defaultParameters
            Raises:
               AssertionError: invalid parameter values supplied
            """
            self._c_param = freesasa_default_parameters
            if param != None:
                  if 'algorithm' in param:    self.setAlgorithm(param['algorithm'])
                  if 'probe-radius' in param: self.setProbeRadius(param['probe-radius'])
                  if 'n-points' in param:     self.setNPoints(param['n-points'])
                  if 'delta' in param:        self.setDelta(param['delta'])
                  if 'n-threads' in param:    self.setNThreads(param['n-threads'])

      def setAlgorithm(self,alg):
            """
            Set algorithm.

            Args:
                alg (str): algorithm name, only allowed values are ShrakeRupley and LeeRichards
            Raises:
                AssertionError: unknwon algorithm specified
            """
            if alg == ShrakeRupley:
                  self._c_param.alg = FREESASA_SHRAKE_RUPLEY
            elif alg == LeeRichards:
                  self._c_param.alg = FREESASA_LEE_RICHARDS
            else:
                  raise AssertionError("Algorithm '%s' is unknown" % alg)

      def algorithm(self):
            """
            Get algorithm.
            
            Returns: 
                Name of algorithm
            """
            if self._c_param.alg == FREESASA_SHRAKE_RUPLEY:
                  return ShrakeRupley
            if self._c_param.alg == FREESASA_LEE_RICHARDS:
                  return LeeRichards
            raise Exception("No algorithm specified, shouldn't be possible")

      def setProbeRadius(self,r):
            """
            Set probe radius.
            
            Args:
                r: probe radius value (>= 0)
            Raises:
                AssertionError: r < 0
            """
            assert(r >= 0)
            self._c_param.probe_radius = r

      def probeRadius(self):
            """
            Get probe radius.

            Returns:
                Probe radius.
            """
            return self._c_param.probe_radius

      def setNPoints(self,n):
            """
            Set number of test points in Shrake & Rupley algorithm.

            Args:
                n (int): Number of points. Must be one of 20, 50, 100, 200, 
                  500, 1000, 2000 or 5000.
            Raises:
                AssertionError: n invalid.
            """
            assert(n in [20,50,100,200,500,1000,2000,5000])
            self._c_param.shrake_rupley_n_points = n

      def nPoints(self):
            """
            Get number of test points in Shrake & Rupley algorithm.

            Returns:
                Number of points.
            """
            return self._c_param.shrake_rupley_n_points

      def setDelta(self,delta):
            """
            Set the value of delta in Lee & Richards algorithm.

            Args:
                delta: Value of delta
            Raises:
                AssertionError: delta must be > 0
            """
            assert(delta > 0)
            self._c_param.lee_richards_delta = delta

      def delta(self):
            """
            Get the value of delta in Lee & Richards algorithm.

            Returns:
                Value of delta.
            """
            return self._c_param.lee_richards_delta

      def setNThreads(self,n):
            """
            Set the number of threads to use in calculations.

            Args:
                n (int): Number of points.
            Raises:
                AssertionError: n muste be > 0.
            """
            assert(n>0)
            self._c_param.n_threads = n

      def nThreads(self):
            """
            Get the number of threads to use in calculations.

            Returns:
                Number of threads.
            """
            return self._c_param.n_threads

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_parameters **p = <freesasa_parameters**> ptr2ptr
            p[0] = &self._c_param

cdef class Result:
      """
      Stores results from SASA calculation.

      The type of object returned by calc(), not intended to be used
      outside of that context. 

      Wraps the C struct freesasa_result.
      """
      cdef freesasa_result* _c_result

      def __cinit__ (self):
            self._c_result = NULL

      def __dealloc__(self):
            if self._c_result is not NULL:
                  freesasa_result_free(self._c_result)

      def nAtoms(self):
            """
            Number of atoms in the results

            Returns:
                Number of atoms
            """
            if self._c_result is not NULL:
                  return self._c_result.n_atoms
            return 0

      def totalArea(self):
            """
            Total SASA.
            
            Returns: 
                The total area in Ångström^2.
            Raises:
                AssertionError: if no results have been associated
                    with the object
            """
            assert(self._c_result is not NULL)
            return self._c_result.total

      def atomArea(self,i):
            """
            SASA for a given atom.

            Args:
                i (int): index of atom.
            Returns:
                SASA of atom i in Ångström^2.
            Raises:
                AssertionError: if no results have been associated
                    with the object or if index is out of bounds
            """
            assert(self._c_result is not NULL)
            assert(i < self._c_result.n_atoms)
            return self._c_result.sasa[i]

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_result **p = <freesasa_result**> ptr2ptr
            p[0] = self._c_result

cdef class Classifier:
      """
      Used to classify atoms by their residue and atom names.

      Wraps the C struct freesasa_classifier. If not initialized with
      filename the default classification is used.

      Residue names should be of the format "ALA","ARG", etc.
      
      Atom names should be of the format " CA ", " N ", etc. The
      default classifier requires these to have length 4, but the
      custom classifiers can handle shorter strings like "CA","N".

      In addition to reading configurations from file, subclasses
      derived from Classifier can be used to define custom atomic
      radii and classes.
      """

      cdef freesasa_classifier* _c_classifier
      
      def __init__ (self,fileName=None):
            """
            Initalizes classifier.
            
            If no file is provided the default classifier is used. The
            config-file should have the same syntax as described in
            the C-documentation and exemplified in share directory.
            
            Args:
                fileName: File with classifier configuration.
            Raises:
                IOError: Problem opening/reading file
                Exception: Problem parsing configuration
            """
            cdef FILE *config
            self._c_classifier = NULL
            if fileName is not None:
                  config = fopen(fileName,'r')
                  if config is NULL:
                        raise IOError("File '%s' could not be opened." % fileName)
                  self._c_classifier = freesasa_classifier_from_file(config)
                  fclose(config)
                  if self._c_classifier is NULL:
                        raise Exception("Error parsing configuration in '%s'." % fileName)

      def __dealloc__(self):
            if self._c_classifier is not NULL:
                  freesasa_classifier_free(self._c_classifier)

      def classify(self,residueName,atomName):
            """
            Class of atom.

            Depending on the configuration these classes can be
            anything, but typically they will be polar/apolar. The
            return values are integers, to turn them into strings, use
            class2str().

            Args:
                residueName (str): Residue name ("ALA","ARG",...).
                atomName (str): Atom name (" CA "," C  ",...).
            Returns:
                An integer representing the class.
            """
            if self._c_classifier is not NULL:
                  return self._c_classifier.sasa_class(residueName,atomName,self._c_classifier)
            return freesasa_default_classifier.sasa_class(residueName,atomName,&freesasa_default_classifier)

      def radius(self,residueName,atomName):
            """
            Radius of atom.

            This allows the classifier to be used to calculate the
            atomic radii used in calculations.

            Args:
                residueName (str): Residue name ("ALA","ARG",...).
                atomName (str): Atom name (" CA "," C  ",...).
            Returns:
                The radius in Ångström^2.
            """
            if self._c_classifier is not NULL:
                  return self._c_classifier.radius(residueName,atomName,self._c_classifier)
            return freesasa_default_classifier.radius(residueName,atomName,&freesasa_default_classifier)

      def class2str(self,classIndex):
            """
            Name of atom class.

            Args:
                classIndex: An integer representing a class, as returned 
                    by classify().
            Returns (str):
                The name of the class
            """
            assert(classIndex >= 0)
            if self._c_classifier is not NULL:
                  assert(classIndex < self._c_classifier.n_classes)
                  return self._c_classifier.class2str(classIndex,self._c_classifier)
            assert(classIndex < freesasa_default_classifier.n_classes)
            return freesasa_default_classifier.class2str(classIndex,&freesasa_default_classifier)

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_classifier **p = <freesasa_classifier**> ptr2ptr
            p[0] = self._c_classifier

# wrapper for freesasa_structure
cdef class Structure:
      """
      Represents a protein structure, including its atomic radii.
      
      Initialized from PDB-file. Calculates atomic radii using default
      classifier, or custom one provided as argument to initalizer

      Wraps the C struct freesasa_structure.
      """

      cdef freesasa_structure* _c_structure
      cdef double* _c_radii

      def __cinit__(self,fileName=None,classifier=None,hetatm=0):
            """
            Initializes Structure

            If PDB file is provided, the structure will be constructed
            based on the file. If not this simply initializes an empty
            structure and the other arguments are ignored. In this case 
            atoms will have to be added manually using addAtom().
            
            Args:
                fileName (str): PDB file (if None empty structure generated)
                classifier: An optional classifier to calculate atomic
                    radii, uses default if none provided
                hetatm (int): Includes hetatms in structure if this is 1
            Raises:
                IOErorr: Problem opening/reading file.
                Exception: Problem parsing PDB file or calculating atomic radii
            """
            self._c_structure = NULL
            self._c_radii = NULL
            if fileName is None:
                  self._c_structure = freesasa_structure_new()
                  return
            cdef FILE *input
            input = fopen(fileName,'r')
            if input is NULL:
                  raise IOError("File '%s' could not be opened." % fileName)
            self._c_structure = freesasa_structure_from_pdb(input,hetatm)
            fclose(input)
            if self._c_structure is NULL:
                  raise Exception("Error reading '%s' as PDB-file." % fileName)
            
            self.setRadiiWithClassifier(classifier)
                        
      def addAtom(self, atomName, residueName, residueNumber, chainLabel, x, y, z):
            """
            Add atom to structure
            
            This function is meant to be used if the structure was
            not initialized from a PDB
            
            Args:
                atomName (str): 4-character string with atom name (e.g. ' CA ')
                residueName (str): 3-character string with residue name (e.g. 'ALA')
                residueNumber (str or int): 4-character string with residue number (e.g. '  12')
                  or integer <= 9999. Some PDBs have residue-numbers that aren't 
                  regular numbers. Therefore treated as a string primarily.
                chainLabel (str): 1-character string with chain label (e.g. 'A')
                x,y,z (float): coordinates
            Raises:
                AssertionError: string-arguments invalid
                Exception: Residue-number invalid
            """
            if (type(residueNumber) is str):
                  resnum = residueNumber
            elif (type(residueNumber) is int):
                  assert residueNumber < 10000
                  resnum = "%4d" % residueNumber
            else:
                  raise Exception("Residue-number invalid, must be either string or number")
            assert len(atomName) == 4
            assert len(residueName) == 3
            assert len(resnum) == 4
            assert len(chainLabel) == 1
            cdef const char *label = chainLabel
            freesasa_structure_add_atom(self._c_structure, atomName,
                                        residueName, resnum, label[0],
                                        x, y, z)

      def setRadiiWithClassifier(self,classifier=None):
            """
            Assign radii to atoms in structure using a classifier

            If no classifier is specified default will be
            used. Over-writes previously assigned radii. Mainly
            intended to be used when atoms were added individually
            using addAtom(). When structure is initialized from PDB a
            classifier can be passed directly to the constructor.
            
            Args: 
                classifier: A classifier to use to calculate radii
            """
            if classifier is None:
                  classifier = Classifier()
            n = self.nAtoms()
            self._c_radii = <double *> realloc(self._c_radii,n*sizeof(double))
            for i in range(0,n):
                  self._c_radii[i] = classifier.radius(self.residueName(i),self.atomName(i))
      
      def setRadii(self,array):
            """
            Set atomic radii from an array

            Args:
                array: Array of atomic radii in Ångström, should 
                    have nAtoms() elements.
            Raises:
                AssertionError: if array has wrong dimension 
            """
            n = self.nAtoms()
            assert len(array) == n
            self._c_radii = <double*> realloc(self._c_radii,n*sizeof(double))
            for i in range(0,n):
                  self._c_radii[i] = array[i]
            
      def nAtoms(self):
            """
            Number of atoms

            Returns:
                Number of atoms
            Raises:
                AssertionError: if not properly initialized
            """
            assert(self._c_structure is not NULL)
            return freesasa_structure_n(self._c_structure)

      def radius(self,i):
            """
            Radius of atom

            Args:
                i (int): Index of atom
            Returns:
                Radius in Ångström
            Raises:
                AssertionError: index out of bounds or object not properly initalized
            """
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_radii is not NULL)
            return self._c_radii[i]
      
      def atomName(self,i):
            """
            Get atom name
            
            Args:
                i (int): Atom index
            Returns:
                Atom name as 4-character string
            Raises:
                AssertionError: if index out of range or 
                   Structure not properly initialized
            """
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            return freesasa_structure_atom_name(self._c_structure,i);
      
      def residueName(self,i):
            """
            Get residue name of given atom
            
            Args:
                i (int): Atom index
            Returns:
                Residue name as 3-character string
            Raises:
                AssertionError: if index out of range or 
                   Structure not properly initialized
            """
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            return freesasa_structure_atom_res_name(self._c_structure,i)

      def residueNumber(self,i):
            """
            Get residue number for given atom
            
            Args:
                i (int): Atom index
            Returns:
                Residue number as 4-character string
            Raises:
                AssertionError: if index out of range or 
                Structure not properly initialized
            """
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            return freesasa_structure_atom_res_number(self._c_structure,i);

      def chainLabel(self,i):
            """
            Get chain label for given atom
            
            Args:
                i (int): Atom index
            Returns:
                Chain label as 1-character string
            Raises:
                AssertionError: if index out of range or 
                Structure not properly initialized
            """
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            cdef char label[2]
            label[0] = freesasa_structure_atom_chain(self._c_structure,i); 
            label[1] = '\0'
            return label

      def _assign_ptr(self, size_t ptr2ptr):
            cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
            p[0] = self._c_structure
      
      def _assign_ptr_radius(self, size_t ptr2ptr):
            cdef double **p = <double**> ptr2ptr
            p[0] = self._c_radii

      def __dealloc__(self):
            if self._c_structure is not NULL:
                  freesasa_structure_free(self._c_structure)
            if self._c_radii is not NULL:
                  free(self._c_radii)

# calculate SASA values for structure

def calc(structure,parameters=None):
      """
      Calculat SASA of Structure

      Args:
          structure: Structure to be used
          parameters: Parameters to use (if not specified defaults are used)
      Returns:
          A Result-object
      Raises:
          Exception: something went wrong in calculation (see C library error messages)
      """
      cdef const freesasa_classifier *c = NULL
      cdef const freesasa_structure *s = NULL
      cdef const freesasa_parameters *p = NULL
      cdef const double *radii = NULL
      if parameters is not None:
            parameters._assign_ptr(<size_t>&p)
      structure._assign_ptr(<size_t>&s)
      structure._assign_ptr_radius(<size_t>&radii)
      r = Result()
      r._c_result = <freesasa_result*> freesasa_calc_structure(s, radii,p)
      if r._c_result is NULL:
            raise Exception("Error calculating SASA.")
      return r
            

# classify results of calculation, returns dictionary with class names
# as keys and areas for the corresponding groups of atoms are values
def classifyResults(result,structure,classifier=None):
      """
      Break SASA result down into classes

      Args:
          result: Result from sasa calculation
          structure: Structure used in calculation
          classifier: Classifier (if not specified default is used)
      Returns:
          Dictionary with name of class and it's SASA value as keys and values
      Raises:
          Exception: Problems with classification, see C library error messages
      """
      if classifier is None:
            classifier = Classifier()
      ret = dict()
      for i in range(0,structure.nAtoms()):
            name = classifier.class2str(classifier.classify(structure.residueName(i),structure.atomName(i))) 
            if name not in ret:
                  ret[name] = 0
            ret[name] += result.atomArea(i)
      return ret

# set verbosity for FreeSASA
def setVerbosity(verbosity):
      """
      Set global verbosity
      
      Args:
          verbosity (int): Can either have values freesasa.silent or freesasa.normal
      Raises:
          AssertionError: if verbosity has illegal value
      """
      assert(verbosity in [silent, normal])
      freesasa_set_verbosity(verbosity)

def getVerbosity():
      """
      Set global verbosity
      
      Returns:
          freesasa.silent or freesasa.normal
      """
      return freesasa_get_verbosity()

