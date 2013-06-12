from ctypes import *
lib = cdll.LoadLibrary('./libpyZipHMM.so')

## HMM IO
def readHMMspec(filename):
    nStates = c_uint()
    nObservables = c_uint()
    lib.c_read_HMM_spec(byref(nStates), byref(nObservables), c_char_p(filename))
    return (nStates, nObservables)

def readHMM(filename):
    pi = Matrix()
    A = Matrix()
    B = Matrix()
    lib.c_read_HMM(pi.obj, A.obj, B.obj, c_char_p(filename))
    return (pi, A, B)

def writeHMM(pi, A, B, filename):
    lib.c_write_HMM(pi.obj, A.obj, B.obj, c_char_p(filename))
    

## Forwarder
lib.Forwarder_new_from_spec.restype = c_void_p
lib.Forwarder_new_from_directory.restype = c_void_p
lib.Forwarder_new_from_directory_and_no_states.restype = c_void_p
lib.Forwarder_forward.restype = c_double
lib.Forwarder_pthread_forward.restype = c_double
lib.Forwarder_get_orig_seq_length.restype = c_uint
lib.Forwarder_get_orig_alphabet_size.restype = c_uint
lib.Forwarder_get_seq_length.restype = c_uint
lib.Forwarder_get_alphabet_size.restype = c_uint
lib.Forwarder_get_pair.restype = py_object


class Forwarder(object):
    def __init__(self, directory = None, nStates = None,
                 seqFilename = None, nObservables = None, nStatesSave = None):
        
        if directory != None and nStates != None:
            self.obj = c_void_p(lib.Forwarder_new_from_directory_and_no_states(c_char_p(directory), nStates))
        elif directory != None:
            self.obj = c_void_p(lib.Forwarder_new_from_directory(c_char_p(directory)))
        elif seqFilename != None and nObservables != None and nStatesSave != None:
            arr = ( c_uint * len(nStatesSave) )()
            arr[:] = nStatesSave
            self.obj = c_void_p(lib.Forwarder_new_from_spec(c_char_p(seqFilename), nObservables, arr, len(nStatesSave)))

    def __del__(self):
        from ctypes import cdll
        lib = cdll.LoadLibrary('./libpyZipHMM.so')
        lib.Forwarder_destructor(self.obj)

    def forward(self, pi, A, B):
        return lib.Forwarder_forward(self.obj, pi.obj, A.obj, B.obj)

    def ptforward(self, pi, A, B, device_filename = None):
        if device_filename == None:
            return lib.Forwarder_pthread_forward(self.obj, pi.obj, A.obj, B.obj, "-")
        else :
            return lib.Forwarder_pthread_forward(self.obj, pi.obj, A.obj, B.obj, device_filename)

    def getOrigSeqLength(self):
        return lib.Forwarder_get_orig_seq_length(self.obj)

    def getOrigAlphabetSize(self):
        return lib.Forwarder_get_orig_alphabet_size(self.obj)

    def getSeqLength(self, no_states):
        return lib.Forwarder_get_seq_length(self.obj, no_states)

    def getAlphabetSize(self, no_states):
        return lib.Forwarder_get_alphabet_size(self.obj, no_states)

    def getPair(self, symbol):
        return lib.Forwarder_get_pair(self.obj, symbol)
         
    def printToDirectory(self, directory):
        lib.Forwarder_print_to_directory(self.obj, c_char_p(directory))

## TimedStopForwarder
lib.TimedStopForwarder_new_from_spec.restype = c_void_p
lib.TimedStopForwarder_new_from_directory.restype = c_void_p
lib.TimedStopForwarder_new_from_directory_and_no_states.restype = c_void_p
lib.TimedStopForwarder_forward.restype = c_double
lib.TimedStopForwarder_pthread_forward.restype = c_double
lib.TimedStopForwarder_get_orig_seq_length.restype = c_uint
lib.TimedStopForwarder_get_orig_alphabet_size.restype = c_uint
lib.TimedStopForwarder_get_seq_length.restype = c_uint
lib.TimedStopForwarder_get_alphabet_size.restype = c_uint
lib.TimedStopForwarder_get_pair.restype = py_object
        
class TimedStopForwarder(object):
    def __init__(self, directory = None, nStates = None,
                 seqFilename = None, nObservables = None, nStatesSave = None, minNoEvals = None):
        
        if directory != None and nStates != None:
            self.obj = c_void_p(lib.TimedStopForwarder_new_from_directory_and_no_states(c_char_p(directory), nStates))
        elif directory != None:
            self.obj = c_void_p(lib.TimedStopForwarder_new_from_directory(c_char_p(directory)))
        elif seqFilename != None and nObservables != None and nStatesSave != None and minNoEvals != None:
            arr = ( c_uint * len(nStatesSave) )()
            arr[:] = nStatesSave
            self.obj = c_void_p(lib.TimedStopForwarder_new_from_spec(c_char_p(seqFilename), nObservables, arr, len(nStatesSave), minNoEvals))
        else:
            raise ValueError("TimedStopForwarder must be given (directory) | (directory, nStates) | (seqFilename, nObservables, nStatesSave, minNoEvals)")

    def __del__(self):
        from ctypes import cdll
        lib = cdll.LoadLibrary('./libpyZipHMM.so')
        lib.TimedStopForwarder_destructor(self.obj)

    def forward(self, pi, A, B):
        return lib.TimedStopForwarder_forward(self.obj, pi.obj, A.obj, B.obj)

    def ptforward(self, pi, A, B, device_filename = None):
        if device_filename == None:
            return lib.TimedStopForwarder_pthread_forward(self.obj, pi.obj, A.obj, B.obj, "-")
        else :
            return lib.TimedStopForwarder_pthread_forward(self.obj, pi.obj, A.obj, B.obj, device_filename)

    def getOrigSeqLength(self):
        return lib.TimedStopForwarder_get_orig_seq_length(self.obj)

    def getOrigAlphabetSize(self):
        return lib.TimedStopForwarder_get_orig_alphabet_size(self.obj)

    def getSeqLength(self, no_states):
        return lib.TimedStopForwarder_get_seq_length(self.obj, no_states)

    def getAlphabetSize(self, no_states):
        return lib.TimedStopForwarder_get_alphabet_size(self.obj, no_states)

    def getPair(self, symbol):
        return lib.TimedStopForwarder_get_pair(self.obj, symbol)
         
    def printToDirectory(self, directory):
        lib.TimedStopForwarder_print_to_directory(self.obj, c_char_p(directory))

## Matrix
lib.Matrix_new_empty.restype = c_void_p
lib.Matrix_new_height_width.restype = c_void_p
lib.Matrix_get_width.restype = c_uint
lib.Matrix_get_height.restype = c_uint
lib.Matrix_get.restype = c_double

class Matrix(object):

    def __init__(self, height = 0, width = 0):
        if height == 0 or width == 0:
            self.obj = c_void_p(lib.Matrix_new_empty())
        else:
            self.obj = c_void_p(lib.Matrix_new_height_width(height, width))

    def __del__(self):
        from ctypes import cdll
        lib = cdll.LoadLibrary('./libpyZipHMM.so')
        lib.Matrix_destructor(self.obj)

    def getWidth(self):
        return lib.Matrix_get_width(self.obj)

    def getHeight(self):
        return lib.Matrix_get_height(self.obj)

    def reset(self, height, width):
        lib.Matrix_reset(self.obj, c_uint(height), c_uint(width))

    def __setitem__(self, (row, column), value):
        lib.Matrix_set(self.obj, c_uint(row), c_uint(column), c_double(value))

    def __getitem__(self, (row, column)):
        return lib.Matrix_get(self.obj, row, column)

    @staticmethod
    def transpose(f, t):
        lib.Matrix_transpose(f.obj, t.obj)

    def p(self):
        lib.Matrix_print(self.obj)

if __name__ == "__main__":
    print "Constructing Matrix(3,7)"
    m = Matrix(3, 7)
    print "Calling getHeight()"
    assert m.getHeight() == 3
    print "Calling getWidth()"
    assert m.getWidth() == 7
    print "Calling setitem method"
    m[1,2] = 0.5
    print "Calling getitem method"
    assert m[1, 2] == 0.5
    print "Calling reset method"
    m.reset(7,3)
    assert m.getHeight() == 7
    assert m.getWidth() == 3

    print "Calling readHMM method"
    (pi, A, B) = readHMM("test_data/test1.hmm")
    assert pi.getHeight() == 2
    assert pi.getWidth()  == 1
    assert A.getHeight()  == 2
    assert A.getWidth()   == 2
    assert B.getHeight()  == 2
    assert B.getWidth()   == 2

    print "Creating Forwarder object from files"
    f = Forwarder(newSeqFilename = "../new_seq.tmp", dataStructureFilename = "../data_structure.tmp")
    assert f.getOrigAlphabetSize() == 2
    assert f.getOrigSeqLength()    == 18
    assert f.getNewAlphabetSize()  == 4
    print "Calling forward on Forwarder object"
    assert abs(f.forward(pi, A, B)  - -12.5671022728) < 0.001

    print "Calling readHMMspec method"
    (nStates, nObservables) = readHMMspec("test_data/test1.hmm")
    assert nStates.value == 2
    assert nObservables.value == 2

    print "Creating Forwarder from sequence and hmm spec"
    f = Forwarder(seqFilename = "test_data/test1.seq", nStates = nStates, nObservables = nObservables)
    assert f.getOrigAlphabetSize() == 2
    assert f.getOrigSeqLength()    == 18
    assert f.getNewAlphabetSize()  == 4
    print "Calling forward"
    assert abs(f.forward(pi, A, B) - -12.5671022728) < 0.001
    
