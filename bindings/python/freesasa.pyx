# Copyright Simon Mitternacht 2013-2015.
#
# This file is part of FreeSASA.
#
# FreeSASA is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# (at your option) any later version.
#
# FreeSASA is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.

## @package freesasa
# @author Simon Mitternacht
# @copyright GPLv3
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
# See documentation of the classes and functions for how to customize behavior
#
# The private methods named _get_address allow access to the C structs
# that these classes are wrapping, which are necessary in the function
# calc(). They are not intended to be used outside of this module.

from libc.stdio cimport FILE, fopen, fclose
from libc.stdlib cimport free, realloc, malloc
from libc.string cimport memcpy
from cpython cimport array
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

## The default values for calculation parameters
defaultParameters = {
     'algorithm'    : ShrakeRupley, 
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
                        raise AssertionError('Key(s): ',unknownKeys,', unknown');

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

      ## The destructor (free internal C objects)
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
#  Wraps the C struct freesasa_classifier. If not initialized with
#  filename the default classification is used.
#
#  Residue names should be of the format `"ALA"`,`"ARG"`, etc.
#      
#  Atom names should be of the format `" CA "`, `" N "`, etc. The
#  default classifier requires these to have length 4, but the
#  custom classifiers can handle shorter strings like `"CA"`,`"N"`.
#
#  In addition to reading configurations from file, subclasses
#  derived from Classifier can be used to define custom atomic
#  radii and/or classes.
cdef class Classifier:
      cdef freesasa_classifier* _c_classifier

      ## Constructor.
      #
      #  If no file is provided the default classifier is used. The
      #  config-file should have the same syntax as described in
      #  the C-documentation and exemplified in share directory.
      #      
      #  @param fileName Name of file with classifier configuration.
      #  @exception IOError   Problem opening/reading file
      #  @exception Exception Problem parsing configuration
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
      
      ## The destructor (free internal C objects)
      def __dealloc__(self):
            if self._c_classifier is not NULL:
                  freesasa_classifier_free(self._c_classifier)

      ## Class of atom.
      #
      #  Depending on the configuration these classes can be
      #  anything, but typically they will be 'polar' and 'apolar'.
      #
      #  @param residueName (str) Residue name (`"ALA"`,`"ARG"`,...).
      #  @param atomName (str) Atom name (`" CA "`,`" C  "`,...).
      #  @return A string describing the class
      def classify(self,residueName,atomName):
            if self._c_classifier is not NULL:
                  classIndex = self._c_classifier.sasa_class(residueName,atomName,self._c_classifier)
                  return self._c_classifier.class2str(classIndex,self._c_classifier)
            classIndex = freesasa_default_classifier.sasa_class(residueName,atomName,self._c_classifier)
            return freesasa_default_classifier.class2str(classIndex,&freesasa_default_classifier)

      ## Radius of atom.
      #
      #  This allows the classifier to be used to calculate the atomic
      #  radii used in calculations.
      #
      #  @param residueName (str) Residue name (`"ALA"`,`"ARG"`,...).
      #  @param atomName (str) Atom name (" CA "," C  ",...).
      #  @return The radius in Å.
      def radius(self,residueName,atomName):
            if self._c_classifier is not NULL:
                  return self._c_classifier.radius(residueName,atomName,self._c_classifier)
            return freesasa_default_classifier.radius(residueName,atomName,&freesasa_default_classifier)

      def _get_address(self, size_t ptr2ptr):
            cdef freesasa_classifier **p = <freesasa_classifier**> ptr2ptr
            p[0] = self._c_classifier

## Represents a protein structure, including its atomic radii.
#      
#  Initialized from PDB-file. Calculates atomic radii using default
#  classifier, or custom one provided as argument to initalizer
#
#  Wraps the C struct freesasa_structure.
cdef class Structure:
      cdef freesasa_structure* _c_structure
      cdef array.array _radii
      ## By default ignore HETATM, Hydrogens and only use first model
      defaultOptions = {'hetatm' : False, 'hydrogen' : False, 'join-models' : False}
      
      ## Constructor
      #
      #  If PDB file is provided, the structure will be constructed
      #  based on the file. If not this simply initializes an empty
      #  structure and the other arguments are ignored. In this case 
      #  atoms will have to be added manually using addAtom().
      #
      #  @param fileName (str) PDB file (if None empty structure generated).
      #  @param classifier An optional Classifier to calculate atomic
      #           radii, uses default if none provided
      #  @param options specify which atoms and models to include
      #  @exception IOError Problem opening/reading file.
      #  @exception Exception Problem parsing PDB file or calculating 
      #                       atomic radii
      def __init__(self,fileName=None,classifier=None,
                   options = defaultOptions):
            self._c_structure = NULL
            self._radii = array.array('d',[])
            if fileName is None:
                  self._c_structure = freesasa_structure_new()
                  return
            cdef FILE *input
            input = fopen(fileName,'r')
            structure_options = Structure._get_structure_options(options)
            if input is NULL:
                  raise IOError("File '%s' could not be opened." % fileName)
            self._c_structure = freesasa_structure_from_pdb(input,structure_options)
            fclose(input)
            if self._c_structure is NULL:
                  raise Exception("Error reading '%s' as PDB-file." % fileName)
            
            self.setRadiiWithClassifier(classifier)
                        
            
      ## Add atom to structure.
      #
      # This function is meant to be used if the structure was
      # not initialized from a PDB
      #      
      # @param atomName (str) 4-character string with atom name (e.g. `" CA "`)
      # @param residueName (str) 3-character string with residue name (e.g. `"ALA"`)
      # @param residueNumber (str or int) 4-character string with residue number (e.g. `'  12'`)
      #      or integer <= 9999. Some PDBs have residue-numbers that aren't 
      #      regular numbers. Therefore treated as a string primarily.
      # @param chainLabel (str) 1-character string with chain label (e.g. 'A')
      # @param x,y,z (float) coordinates
      # 
      # @exception AssertionError string-arguments invalid
      # @exception Exception Residue-number invalid
      def addAtom(self, atomName, residueName, residueNumber, chainLabel, x, y, z):
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

      ## Assign radii to atoms in structure using a classifier.
      #
      # If no classifier is specified default will be
      # used. Over-writes previously assigned radii. Mainly intended
      # to be used when atoms were added individually using
      # addAtom(). When structure is initialized from PDB a classifier
      # can be passed directly to the constructor.
      #
      # @param classifier A classifier to use to calculate radii
      def setRadiiWithClassifier(self,classifier=None):
            if classifier is None:
                  classifier = Classifier()
            n = self.nAtoms()
            array.resize(self._radii,n*sizeof(double))
            for i in range(0,n):
                  self._radii[i] = classifier.radius(self.residueName(i),self.atomName(i))
      
      ## Set atomic radii from an array
      # @param radiusArray: Array of atomic radii in Ångström, should 
      #                     have nAtoms() elements.
      # @exception AssertionError if radiusArray has wrong dimension 
      def setRadii(self,radiusArray):
            n = self.nAtoms()
            assert len(radiusArray) == n
            array.resize(self._radii,n*sizeof(double))
            for i in range(0,n):
                  self._radii[i] = radiusArray[i]
            
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
      # @exception AssertionError if index out of bounds or object not properly initalized
      def radius(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._radii.ob_size > 0)
            return self._radii[i]
      
      ## Get atom name
      # @param i (int) Atom index.
      # @return Atom name as 4-character string.
      # @exception AssertionError: if index out of range or Structure not properly initialized.
      def atomName(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            return freesasa_structure_atom_name(self._c_structure,i);
      
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
            return freesasa_structure_atom_res_number(self._c_structure,i);

      ## Get chain label for given atom.
      # @param i (int) Atom index.
      # @return Chain label as 1-character string
      # @exception AssertionError if index out of range or Structure not properly initialized
      def chainLabel(self,i):
            assert(i >= 0 and i < self.nAtoms())
            assert(self._c_structure is not NULL)
            cdef char label[2]
            label[0] = freesasa_structure_atom_chain(self._c_structure,i); 
            label[1] = '\0'
            return label

      @staticmethod
      def _get_structure_options(param):
            options = 0
            
            # check validity of options
            knownOptions = {'hetatm','hydrogen','join-models','separate-models','separate-chains'}
            unknownOptions = []
            for key in param:
                  if not key in knownOptions:
                        unknownOptions.append(key)
            if len(unknownOptions) > 0:
                  raise AssertionError("Option(s): ",unknownOptions," unknown.");

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
            return options

      def _get_address(self, size_t ptr2ptr):
            cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
            p[0] = self._c_structure

      def _get_address_radius(self, size_t ptr2ptr):
            cdef double **p = <double**> ptr2ptr
            p[0] = self._radii.data.as_doubles

      def _set_address(self, size_t ptr2ptr):
            cdef freesasa_structure **p = <freesasa_structure**> ptr2ptr
            self._c_structure = p[0]

      ## The destructor (free internal C objects)
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
# @exception IOError if can't open file
# @exception Exception: if there are problems parsing the input
# @return An array of Structures
def structureArray(fileName,
                   options = defaultStructureArrayOptions,
                   classifier = None):
      assert fileName is not None
      structure_options = Structure._get_structure_options(options)
      cdef FILE *input
      input = fopen(fileName,'r')
      if input is NULL:
            raise IOError("File '%s' could not be opened." % fileName)
      cdef int n;
      cdef freesasa_structure** sArray = freesasa_structure_array(input,&n,structure_options)
      fclose(input)
      if sArray is NULL:
            raise Exception("Problems reading structures in '%s'." % fileName)
      structures = []
      for i in range(0,n):
            structures.append(Structure())
            structures[-1]._set_address(<size_t> &sArray[i])
            structures[-1].setRadiiWithClassifier(classifier)
      free(sArray)
      return structures


## Calculate SASA of Structure
# @param structure Structure to be used
# @param parameters Parameters to use (if not specified defaults are used)
# @return A Result object
# @exception Exception something went wrong in calculation (see C library error messages)
def calc(structure,parameters=None):
      cdef const freesasa_parameters *p = NULL
      cdef const freesasa_structure *s = NULL
      cdef const double *radii = NULL
      if parameters is not None:  parameters._get_address(<size_t>&p)
      structure._get_address(<size_t>&s)
      structure._get_address_radius(<size_t>&radii)
      result = Result()
      result._c_result = <freesasa_result*> freesasa_calc_structure(s, radii,p)
      if result._c_result is NULL:
            raise Exception("Error calculating SASA.")
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
#   syntax, e.g. "s1, resn ala+arg" or "s2, chain A and resi 1-5".
# @param structure A Structure.  
# @param result Result from sasa calculation on structure.
# @return Dictionary with names of selections ("s1","s2",...) as 
#   keys, and the corresponding SASA values as values.
# @exception Exception: Parser failed (typically syntax error), see
#   library error messages.
def selectArea(commands,structure,result):
      cdef freesasa_structure *s
      cdef freesasa_result *r
      cdef double area
      cdef char *name = <char*>malloc(FREESASA_MAX_SELECTION_NAME+1);
      structure._get_address(<size_t> &s)
      result._get_address(<size_t> &r)
      value = dict()
      for cmd in commands:
            ret = freesasa_select_area(cmd,name,&area,s,r)
            if ret == FREESASA_FAIL:
                  raise Exception("Error parsing '%s'" % cmd)
            value[name] = area
      free(name)
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

