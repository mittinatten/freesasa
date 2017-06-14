## @package freesasa
# @author Simon Mitternacht
# @copyright [MIT License](md_license.html)
#
# @brief FreeSASA Python interface.
#
# A minimal program would be something like
#
# ~~~{.py}
#     structure = freesasa.Structure("1abc.pdb")
#     result = freesasa.calc(structure)
#     print "SASA = %f" % result.totalArea()
# ~~~
#
# See documentation of the classes and functions for how to customize behavior.
#

# cython: c_string_type=str, c_string_encoding=ascii

from libc.stdio cimport FILE, fopen, fclose
from libc.stdlib cimport free, realloc, malloc
from libc.string cimport memcpy
from cfreesasa cimport *

## Used to specify the algorithm by Shrake & Rupley
ShrakeRupley = 'ShrakeRupley'

## Used to specify the algorithm by Lee & Richards
LeeRichards = 'LeeRichards'

## Used for classification
polar = 'Polar'

## Used for classification
apolar = 'Apolar'

## int: Suppress all warnings and errors (used by setVerbosity())
silent = FREESASA_V_SILENT

## int: Suppress all warnings but not errors (used by setVerbosity())
nowarnings = FREESASA_V_NOWARNINGS

## int: Normal verbosity (used by setVerbosity())
normal = FREESASA_V_NORMAL

## int: Print debug messages (used by setVerbosity())
debug = FREESASA_V_DEBUG

## The default values for calculation parameters
defaultParameters = {
     'algorithm'    : LeeRichards, 
     'probe-radius' : freesasa_default_parameters.probe_radius,
     'n-points'     : freesasa_default_parameters.shrake_rupley_n_points,
     'n-slices'     : freesasa_default_parameters.lee_richards_n_slices,
     'n-threads'    : freesasa_default_parameters.n_threads 
}

## Stores parameter values to be used by calculation.
#   
#  Wraps the C struct freesasa_parameters.
cdef class Parameters:

      cdef freesasa_parameters _c_param
      ## Initializes Parameters object.
      #      
      #  @param param (dict) optional argument to specify parameter-values, modeled after
      #            ::defaultParameters
      #  @exception AssertionError invalid parameter values supplied
      #  @see defaultParameters
      def __init__(self,param=None):
   
            self._c_param = freesasa_default_parameters
            if param != None:
                  if 'algorithm' in param:    self.setAlgorithm(param['algorithm'])
                  if 'probe-radius' in param: self.setProbeRadius(param['probe-radius'])
                  if 'n-points' in param:     self.setNPoints(param['n-points'])
                  if 'n-slices' in param:     self.setNSlices(param['n-slices'])
                  if 'n-threads' in param:    self.setNThreads(param['n-threads'])
                  unknownKeys = []
                  for key in param:
                        if not key in defaultParameters:
                              unknownKeys.append(key)
                  if len(unknownKeys) > 0:
                        raise AssertionError('Key(s): ',unknownKeys,', unknown')

      ## Set algorithm.
      #
      #  @param alg (str) algorithm name, only allowed values are ::ShrakeRupley and ::LeeRichards
      #  @exception AssertionError unknown algorithm specified
      def setAlgorithm(self,alg):
            if alg == ShrakeRupley:
                  self._c_param.alg = FREESASA_SHRAKE_RUPLEY
            elif alg == LeeRichards:
                  self._c_param.alg = FREESASA_LEE_RICHARDS
            else:
                  raise AssertionError("Algorithm '%s' is unknown" % alg)

      ## Get algorithm.
      #     
      #  @return Name of algorithm
      def algorithm(self):
            if self._c_param.alg == FREESASA_SHRAKE_RUPLEY:
                  return ShrakeRupley
            if self._c_param.alg == FREESASA_LEE_RICHARDS:
                  return LeeRichards
            raise Exception("No algorithm specified, shouldn't be possible")

      ## Set probe radius.
      # @param r probe radius in Å (>= 0)
      # @exception AssertionError r < 0
      def setProbeRadius(self,r):
            assert(r >= 0)
            self._c_param.probe_radius = r

      ## Get probe radius.
      #  @return Probe radius in Å
      def probeRadius(self):
            return self._c_param.probe_radius

      ## Set number of test points in Shrake & Rupley algorithm.
      #  @param n (int) Number of points (> 0). 
      #  @exception AssertionError n <= 0.
      def setNPoints(self,n):
            assert(n > 0)
            self._c_param.shrake_rupley_n_points = n

      ## Get number of test points in Shrake & Rupley algorithm.
      #  @return Number of points.
      def nPoints(self):
            return self._c_param.shrake_rupley_n_points

      ## Set the number of slices per atom in Lee & Richards algorithm.
      #  @param n (int) Number of slices (> 0)
      #  @exception AssertionError n <= 0
      def setNSlices(self,n):
            assert(n> 0)
            self._c_param.lee_richards_n_slices = n

      ## Get the number of slices per atom in Lee & Richards algorithm.
      #  @return Number of slices.
      def nSlices(self):
            return self._c_param.lee_richards_n_slices

      
      ## Set the number of threads to use in calculations.
      #  @param n (int) Number of points (> 0)
      #  @exception AssertionError n <= 0
      def setNThreads(self,n):
            assert(n>0)
            self._c_param.n_threads = n

            
      ## Get the number of threads to use in calculations.
      #  @return Number of threads.
      def nThreads(self):
            return self._c_param.n_threads

      # not pretty, but only way I've found to pass pointers around
      def _get_address(self, size_t ptr2ptr):
            cdef freesasa_parameters **p = <freesasa_parameters**> ptr2ptr
            p[0] = &self._c_param

## Stores results from SASA calculation.
#  The type of object returned by calc(), not intended to be used
#  outside of that context. 
#
#  Wraps the C struct freesasa_result.
cdef class Result:
      cdef freesasa_result* _c_result

      ## The constructor
      def __cinit__ (self):
            self._c_result = NULL

      ## The destructor
      def __dealloc__(self):
            if self._c_result is not NULL:
                  freesasa_result_free(self._c_result)
            
      ## Number of atoms in the results
      #  @return Number of atoms
      def nAtoms(self):
            if self._c_result is not NULL:
                  return self._c_result.n_atoms
            return 0

      
      ## Total SASA.
      # @return The total area in Å^2.
      # @exception AssertionError If no results have been associated
      #              with the object
      def totalArea(self):
            assert(self._c_result is not NULL)
            return self._c_result.total

      ## SASA for a given atom.
      #  @param i (int) index of atom.
      #  @return SASA of atom i in Å^2.
      #  @exception AssertionError If no results have been associated
      #              with the object or if index is out of bounds
      def atomArea(self,i):
            assert(self._c_result is not NULL)
            assert(i < self._c_result.n_atoms)
            return self._c_result.sasa[i]

      def _get_address(self, size_t ptr2ptr):
            cdef freesasa_result **p = <freesasa_result**> ptr2ptr
            p[0] = self._c_result

      
## Assigns class and radius to atom by residue and atom name.
#
#  Subclasses derived from Classifier can be used to define custom
#  atomic radii and/or classes. Can also be initialized from a @ref
#  Config-file "configuration file" with a custom classifier.
#
#  Wraps a C freesasa_classifier. If initialized without arguments the
#  default classifier is used.
#
#  Derived classifiers should set the member purePython to True
#
#  Residue names should be of the format `"ALA"`,`"ARG"`, etc.
#      
#  Atom names should be of the format `"CA"`, `"N"`, etc. 
#
cdef class Classifier:
      cdef freesasa_classifier* _c_classifier
      purePython = False

      ## Constructor.
      #
      #  If no file is provided the default classifier is used. 
      #
      #  @see @ref Config-file.
      #      
      #  @param fileName Name of file with classifier configuration.
      #  @exception IOError   Problem opening/reading file
      #  @exception Exception Problem parsing provided configuration or 
      #                       initializing defaults
      def __cinit__ (self,fileName=None):
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
            else:
                  self._c_classifier = &freesasa_default_classifier
      ## The destructor
      def __dealloc__(self):
            if self._c_classifier is not &freesasa_default_classifier:
                  freesasa_classifier_free(self._c_classifier)


      # This is used internally to determine if a Classifier wraps a C
      # classifier or not (necessary when generating structures)
      # @return Boolean
      def _isCClassifier(self):
            return not self.purePython

      ## Class of atom.
      #
      #  Depending on the configuration these classes can be
      #  anything, but typically they will be 'Polar' and 'Apolar'. 
      #  Unrecognized atoms will get the class 'Unknown'. 
      #
      #  @param residueName (str) Residue name (`"ALA"`,`"ARG"`,...).
      #  @param atomName (str) Atom name (`"CA"`,`"C"`,...).
      #  @return A string describing the class
      def classify(self,residueName,atomName):
            classIndex = freesasa_classifier_class(self._c_classifier, residueName, atomName)
            return freesasa_classifier_class2str(classIndex)

      ## Radius of atom.
      #
      #  This allows the classifier to be used to calculate the atomic
      #  radii used in calculations. Unknown atoms will get a negative
      #  radius.
      #
      #  @param residueName (str) Residue name (`"ALA"`, `"ARG"`, ...).
      #  @param atomName (str) Atom name (`"CA"`, `"C"`, ...).
      #  @return The radius in Å.
      def radius(self,residueName,atomName):
            return freesasa_classifier_radius(self._c_classifier, residueName, atomName)

      def _get_address(self, size_t ptr2ptr):
            cdef freesasa_classifier **p = <freesasa_classifier**> ptr2ptr
            p[0] = self._c_classifier

## Represents a protein structure, including its atomic radii.
#      
#  Initialized from PDB-file. Calculates atomic radii using default
#  classifier, or custom one provided as argument to initalizer
#
#  Wraps the C struct freesasa_structure.
#
#  Since it is intended to be a static structure the word 'get' is
#  omitted in the getter-functions.
cdef class Structure:
      cdef freesasa_structure* _c_structure
      ## By default ignore HETATM, Hydrogens, only use first model. For unknown atoms
      # try to guess the radius, if this fails, assign radius 0 (to
      # allow changing the radius later).
      defaultOptions = {         'hetatm' : False,
            'hydrogen' : False,
            'join-models' : False,
            'skip-unknown' : False,
            'halt-at-unknown' : False
            }

      ## Constructor
      #
      #  If PDB file is provided, the structure will be constructed
      #  based on the file. If not, this simply initializes an empty
      #  structure and the other arguments are ignored. In this case 
      #  atoms will have to be added manually using addAtom().
      #
      #  @param fileName (str) PDB file (if None empty structure generated).
      #  @param classifier An optional Classifier to calculate atomic
      #           radii, uses default if none provided
      #  @param options specify which atoms and models to include
      #  @exception IOError Problem opening/reading file.
      #  @exception Exception Problem parsing PDB file or calculating 
      #      atomic radii.
      #  @exception Exception If option 'halt-at-unknown' selected and
      #      unknown atom encountered.
      def __init__(self,fileName=None,classifier=None,
                   options = defaultOptions):

            self._c_structure = NULL
            cdef freesasa_classifier *c = NULL
            if classifier is None:
                  classifier = Classifier()
            if classifier._isCClassifier():
                  classifier._get_address(<size_t>&c)

            if fileName is None:
                  self._c_structure = freesasa_structure_new()
                  return
            cdef FILE *input
            input = fopen(fileName,'r')
            if input is NULL:
                  raise IOError("File '%s' could not be opened." % fileName)
            structure_options = Structure._get_structure_options(options)

            if not classifier._isCClassifier(): # supress warnings
                  setVerbosity(silent)

            self._c_structure = freesasa_structure_from_pdb(input, c, structure_options)

            if not classifier._isCClassifier():
                  setVerbosity(normal)

            fclose(input)

            if self._c_structure is NULL:
                  raise Exception("Error reading '%s'." % fileName)

            # for pure Python classifiers we use the default
            # classifier above to initialize the structure and then
            # reassign radii using the provided classifier here
            if (not classifier._isCClassifier()):
                  self.setRadiiWithClassifier(classifier)


      ## Add atom to structure.
      #
      # This function is meant to be used if the structure was not
      # initialized from a PDB. Default radii will be assigned to each
      # atom. This can be overriden by calling
      # setRadiiWithClassifier() afterwards.
      #
      # There are no restraints on string lengths for the arguments, but
      # the atom won't be added if the default classifier doesn't
      # recognize the atom and also cannot deduce its element from the
      # atom name.
      #      
      # @param atomName (str) atom name (e.g. `"CA"`)
      # @param residueName (str) residue name (e.g. `"ALA"`)
      # @param residueNumber (str or int) residue number (e.g. `'12'`)
      #      or integer. Some PDBs have residue-numbers that aren't 
      #      regular numbers. Therefore treated as a string primarily.
      # @param chainLabel (str) 1-character string with chain label (e.g. 'A')
      # @param x,y,z (float) coordinates
      # 
      # @exception Exception Residue-number invalid
      def addAtom(self, atomName, residueName, residueNumber, chainLabel, x, y, z):
            if (type(residueNumber) is str):
                  resnum = residueNumber
            elif (type(residueNumber) is int):
                  resnum = "%d" % residueNumber
            else:
                  raise Exception("Residue-number invalid, must be either string or number")
            cdef const char *label = chainLabel
            ret = freesasa_structure_add_atom(self._c_structure, atomName,
                                              residueName, resnum, label[0],
                                              x, y, z)
            assert(ret != FREESASA_FAIL)

      ## Assign radii to atoms in structure using a classifier.
      #
      # @param classifier A classifier to use to calculate radii
      # @exception AssertionError if structure not properly initialized
      def setRadiiWithClassifier(self,classifier):
            assert(self._c_structure is not NULL)
            n = self.nAtoms()
            r = []
            for i in range(0,n):
                  r.append(classifier.radius(self.residueName(i), self.atomName(i)))
            self.setRadii(r)

      ## Set atomic radii from an array
      # @param radiusArray Array of atomic radii in Ångström, should 
      #                    have nAtoms() elements.
      # @exception AssertionError if radiusArray has wrong dimension, structure 
      #                           not properly initialized, or if the array contains
      #                           negative radii (not properly classified?) 
      def setRadii(self,radiusArray):
            assert(self._c_structure is not NULL)
            n = self.nAtoms()
            assert len(radiusArray) == n
            cdef double *r = <double *>malloc(sizeof(double)*n)
            assert(r is not NULL)
            for i in range(0,n):
                  r[i] = radiusArray[i]
                  assert(r[i] >= 0)
            freesasa_structure_set_radius(self._c_structure, r)

      ## Number of atoms.
      #
      # @return  Number of atoms
      # @exception AssertionError if not properly initialized
      def nAtoms(self):
            assert(self._c_structure is not NULL)
            return freesasa_structure_n(self._c_structure)

      ## Radius of atom.
      # @param i (int) Index of atom.
      # @return Radius in Å. 
      # @exception AssertionError if index out of bounds, object not properly initalized.
      def radius(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            cdef const double *r = freesasa_structure_radius(self._c_structure)
            assert(r is not NULL)
            return r[i]

      ## Set radius for a given atom
      # @param atomIndex Index of atom
      # @param radius Value of radius
      # @exception AssertionError if index out of bounds, radius
      #            negative, or structure not properly initialized
      def setRadius(self, atomIndex, radius):
            assert(self._c_structure is not NULL)
            assert(atomIndex >= 0 and atomIndex < self.nAtoms())
            assert(radius >= 0)
            freesasa_structure_atom_set_radius(self._c_structure, atomIndex, radius)
      
      ## Get atom name
      # @param i (int) Atom index.
      # @return Atom name as 4-character string.
      # @exception AssertionError: if index out of range or Structure not properly initialized.
      def atomName(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            return freesasa_structure_atom_name(self._c_structure,i)
      
      ## Get residue name of given atom.
      # @param i (int) Atom index.
      # @return Residue name as 3-character string.
      # @exception AssertionError if index out of range or Structure not properly initialized
      def residueName(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            return freesasa_structure_atom_res_name(self._c_structure,i)

      ## Get residue number for given atom
      # @param i (int) Atom index.
      # @return Residue number as 4-character string
      # @exception AssertionError if index out of range or Structure not properly initialized
      def residueNumber(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            return freesasa_structure_atom_res_number(self._c_structure,i)

      ## Get chain label for given atom.
      # @param i (int) Atom index.
      # @return Chain label as 1-character string
      # @exception AssertionError if index out of range or Structure not properly initialized
      def chainLabel(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            cdef char label[2]
            label[0] = freesasa_structure_atom_chain(self._c_structure,i) 
            label[1] = '\0'
            return label

      ## Get coordinates of given atom.
      # @param i (int) Atom index.
      # @return array of x, y, and z coordinates
      # @exception AssertionError if index out of range or Structure not properly initialized
      def coord(self, i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            cdef const double *coord = freesasa_structure_coord_array(self._c_structure)
            return [coord[3*i], coord[3*i+1], coord[3*i+2]]

      @staticmethod
      def _get_structure_options(param):
            options = 0
            
            # check validity of options
            knownOptions = {'hetatm','hydrogen','join-models','separate-models',
                            'separate-chains','skip-unknown','halt-at-unknown'}
            unknownOptions = []
            for key in param:
                  if not key in knownOptions:
                        unknownOptions.append(key)
            if len(unknownOptions) > 0:
                  raise AssertionError("Option(s): ",unknownOptions," unknown.")

            # calculate bitfield
            if 'hetatm' in param and param['hetatm']:
                  options |= FREESASA_INCLUDE_HETATM
            if 'hydrogen' in param and param['hydrogen']:
                  options |= FREESASA_INCLUDE_HYDROGEN
            if 'join-models' in param and param['join-models']:
                  options |= FREESASA_JOIN_MODELS
            if 'separate-models' in param and param['separate-models']:
                  options |= FREESASA_SEPARATE_MODELS
            if 'separate-chains' in param and param['separate-chains']:
                  options |= FREESASA_SEPARATE_CHAINS
            if 'skip-unknown' in param and param['skip-unknown']:
                  options |= FREESASA_SKIP_UNKNOWN
            if 'halt-at-unknown' in param and param['halt-at-unknown']:
                  options |= FREESASA_HALT_AT_UNKNOWN
            return options

      def _get_address(self, size_t ptr2ptr):
            cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
            p[0] = self._c_structure

      def _set_address(self, size_t ptr2ptr):
            cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
            self._c_structure = p[0]

      ## The destructor
      def __dealloc__(self):
            if self._c_structure is not NULL:
                  freesasa_structure_free(self._c_structure)

## Default options for structureArray.
# Defined separately for Doxygen's sake.
defaultStructureArrayOptions = {
      'hetatm' : False, 
      'hydrogen' : False,
      'separate-chains' : True,
      'separate-models' : False
}

## Create array of structures from PDB file.
#
# Split PDB file into several structures by either by treating
# chains separately, by treating each MODEL as a separate
# structure, or both.
# 
# @param fileName (str) The PDB file.
# @param options (dict) Specification for how to read the PDB-file
#  (see def value for options).
# @param classifier: Classifier to assign atoms radii, default is used
#   if none specified.
# @exception AssertionError if fileName is None
# @exception AssertionError if an option value is not recognized
# @exception AssertionError if neither of the options 'separate-chains'
#   and 'separate-models' are specified.
# @exception IOError if can't open file
# @exception Exception: if there are problems parsing the input
# @return An array of Structures
def structureArray(fileName,
                   options = defaultStructureArrayOptions,
                   classifier = None):
      assert fileName is not None
      # we need to have at least one of these
      assert(('separate-chains' in options and options['separate-chains'] is True) 
             or ('separate-models' in options and options['separate-models'] is True))
      structure_options = Structure._get_structure_options(options)
      cdef FILE *input
      input = fopen(fileName,'r')
      if input is NULL:
            raise IOError("File '%s' could not be opened." % fileName)
      cdef int n
      cdef freesasa_structure** sArray = freesasa_structure_array(input,&n,NULL,structure_options)
      fclose(input)
      if sArray is NULL:
            raise Exception("Problems reading structures in '%s'." % fileName)
      structures = []
      for i in range(0,n):
            structures.append(Structure())
            structures[-1]._set_address(<size_t> &sArray[i])
            if classifier is not None:
                  structures[-1].setRadiiWithClassifier(classifier)
      free(sArray)
      return structures


## Calculate SASA of Structure
# @param structure Structure to be used
# @param parameters Parameters to use (if not specified defaults are used)
# @return A Result object
# @exception Exception: something went wrong in calculation (see C library error messages)
def calc(structure,parameters=None):
      cdef const freesasa_parameters *p = NULL
      cdef const freesasa_structure *s = NULL
      if parameters is not None:  parameters._get_address(<size_t>&p)
      structure._get_address(<size_t>&s)
      result = Result()
      result._c_result = <freesasa_result*> freesasa_calc_structure(s,p)
      if result._c_result is NULL:
            raise Exception("Error calculating SASA.")
      return result

## Calculate SASA for a set of coordinates and radii
# @param coord list of size 3*N with atomic coordinates (x1, y1, z1,
#   x2, y2, z2, ..., x_N, y_N, z_N'.
# @param radii array of size N with atomic radii (r_1, r_2, ..., r_N)
# @param Parameters to use (if not specified, defaults are used)
# @exception AssertionError: mismatched array-sizes
# @exception Exception: Out of memory
# @exception Exception: something went wrong in calculation (see C library error messages)
def calcCoord(coord, radii, parameters=None):
      assert(len(coord) == 3*len(radii))

      cdef const freesasa_parameters *p = NULL
      cdef double *c = <double*> malloc(len(coord)*sizeof(double))
      cdef double *r = <double*> malloc(len(radii)*sizeof(double))
      if c is NULL or r is NULL:
            raise Exception("Memory allocation error")

      for i in xrange(len(coord)):
            c[i] = coord[i]
      for i in xrange(len(radii)):
            r[i] = radii[i]

      if parameters is not None: parameters._get_address(<size_t>&p)

      result = Result()
      result._c_result = <freesasa_result*> freesasa_calc_coord(c, r, len(radii), p)

      if result._c_result is NULL:
            raise Exception("Error calculating SASA.")

      free(c)
      free(r)

      return result

## Break SASA result down into classes.
# @param result Result from sasa calculation.
# @param structure Structure used in calculation.
# @param classifier Classifier (if not specified default is used).
# @return Dictionary with names of classes as keys and their SASA values as values.
# @exception Exception: Problems with classification, see C library error messages 
#  (or Python exceptions if run with derived classifier).
def classifyResults(result,structure,classifier=None):
      if classifier is None:
            classifier = Classifier()
      ret = dict()
      for i in range(0,structure.nAtoms()):
            name = classifier.classify(structure.residueName(i),structure.atomName(i))
            if name not in ret:
                  ret[name] = 0
            ret[name] += result.atomArea(i)
      return ret

## Sum SASA result over a selection of atoms
# @param commands A list of commands with selections using Pymol
#   syntax, e.g. `"s1, resn ala+arg"` or `"s2, chain A and resi 1-5"` 
#   (see @ref Selection).
# @param structure A Structure.  
# @param result Result from sasa calculation on structure.
# @return Dictionary with names of selections (`"s1"`,`"s2"`,...) as 
#   keys, and the corresponding SASA values as values.
# @exception Exception: Parser failed (typically syntax error), see
#   library error messages.
def selectArea(commands, structure, result):
      cdef freesasa_structure *s
      cdef freesasa_result *r
      cdef freesasa_selection *selection
      structure._get_address(<size_t> &s)
      result._get_address(<size_t> &r)
      value = dict()
      for cmd in commands:
            selection = freesasa_selection_new(cmd, s, r)
            if selection == NULL:
                  raise Exception("Error parsing '%s'" % cmd)
            value[freesasa_selection_name(selection)] = freesasa_selection_area(selection)
            freesasa_selection_free(selection)
      return value

## Set global verbosity
# @param verbosity Can have values freesasa.silent, freesasa.nowarnings or freesasa.normal
# @exception AssertionError if verbosity has illegal value
def setVerbosity(verbosity):
      assert(verbosity in [silent, nowarnings, normal])
      freesasa_set_verbosity(verbosity)


## Get global verbosity
# @return Verbosity (freesasa.silent, freesasa.nowarnings or freesasa.normal)
def getVerbosity():
      return freesasa_get_verbosity()

## Create a freesasa structure from a Bio.PDB structure
#
#  @remark Experimental, not thorougly tested yet
#
#  @param bioPDBStructure a Bio.PDB structure
#  @param classifier an optional classifier to specify atomic radii
#  @param options Options supported are 'hetatm', 'skip-unknown' and 'halt-at-unknown'
#  @return a freesasa.Structure
#
#  @exception Exception if option 'halt-at-unknown' is selected and
#             unknown atoms are encountered. Passes on exceptions from
#             Structure.addAtom() and
#             Structure.setRadiiWithClassifier().
def structureFromBioPDB(bioPDBStructure, classifier=None, options = Structure.defaultOptions):
      structure = Structure()
      if (classifier is None):
            classifier = Classifier()
      optbitfield = Structure._get_structure_options(options)

      atoms = bioPDBStructure.get_atoms()

      for a in atoms:
            r = a.get_parent()
            hetflag, resseq, icode = r.get_id()

            if (hetflag is not ' ' and not (optbitfield & FREESASA_INCLUDE_HETATM)):
                  continue

            c = r.get_parent()
            v = a.get_vector()

            if (classifier.classify(r.get_resname(), a.get_fullname()) is 'Unknown'):
                  if (optbitfield & FREESASA_SKIP_UNKNOWN):
                        continue
                  if (optbitfield & FREESASA_HALT_AT_UNKNOWN):
                        raise Exception("Halting at unknown atom")

            structure.addAtom(a.get_fullname(), r.get_resname(), resseq, c.get_id(),
                              v[0], v[1], v[2])

      structure.setRadiiWithClassifier(classifier)
      return structure

## Calc SASA from Bio.PDB structure
#
#  Usage 
#
#      result, sasa_classes = calcBioPDB(structure, ...)  
#
#  @remark Experimental, not thorougly tested yet
#
#  @param bioPDBStructure A Bio.PDB structure
#  @param parameters A freesasa.Paremeters object
#  @param classifier A freesasa.Classifier object
#  @param options Options supported are 'hetatm', 'skip-unknown' and 'halt-at-unknown'
#
#  @return A freesasa.Result object and a dictionary with classes
#         defined by the classifier and associated areas
#
#  @exception Exception if unknown atom is encountered and the option
#             'halt-at-unknown' is active. Passes on exceptions from
#             calc(), classifyResults() and structureFromBioPDB().
def calcBioPDB(bioPDBStructure, parameters = Parameters(), 
               classifier = None, options = Structure.defaultOptions):
      structure = structureFromBioPDB(bioPDBStructure, classifier, options)
      result = calc(structure, parameters)
      sasa_classes = classifyResults(result, structure, classifier)
      return result, sasa_classes
      

