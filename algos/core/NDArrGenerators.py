'''
Created on Nov 23, 2015

@author: Mikhail
'''
import numpy as np
#-----------------------------

def prod_of_elements(arr, dtype=np.int) :
    """Returns product of sequence elements
    """
    return np.prod(arr,axis=None,dtype=dtype)
#-----------------------------

def size_from_shape(shape) :
    """Returns size from the shape sequence 
    """
    return prod_of_elements(shape)

#-----------------------------

def shape_as_2d(sh) :
    """Returns 2-d shape for n-d shape if n>2, otherwise returns unchanged shape.
    """
    if len(sh)<3 : return sh
    return (size_from_shape(sh)/sh[-1], sh[-1])

#-----------------------------

def shape_as_3d(sh) :
    """Returns 3-d shape for n-d shape if n>3, otherwise returns unchanged shape.
    """
    if len(sh)<4 : return sh
    return (size_from_shape(sh)/sh[-1]/sh[-2], sh[-2], sh[-1])

#-----------------------------

def reshape_to_2d(arr) :
    """Returns n-d re-shaped to 2-d
    """
    arr.shape = shape_as_2d(arr.shape)
    return arr

#-----------------------------

def reshape_to_3d(arr) :
    """Returns n-d re-shaped to 3-d
    """
    arr.shape = shape_as_3d(arr.shape)
    return arr

#-----------------------------

def random_standard(shape=(40,60), mu=200, sigma=25) :
    """Returns numpy array of requested shape and type filled with normal distribution for mu and sigma.
    """
    return mu + sigma*np.random.standard_normal(shape)

#-----------------------------

def random_exponential(shape=(40,60), a0=100) :
    """Returns numpy array of requested shape and type filled with exponential distribution for width a0.
    """
    return a0*np.random.standard_exponential(size=shape)

#-----------------------------

def random_1(shape=(40,60), dtype=np.float) :
    """Returns numpy array of requested shape and type filled with random numbers in the range [0,255].
    """
    a = np.random.random(shape)
    return np.require(a, dtype) 

#-----------------------------

def random_256(shape=(40,60), dtype=np.uint8) :
    """Returns numpy array of requested shape and type filled with random numbers in the range [0,255].
    """
    a = 255*np.random.random(shape)
    return np.require(a, dtype) 

#-----------------------------

def random_xffffffff(shape=(40,60), dtype=np.uint32, add=0xff000000) :
    """Returns numpy array of requested shape and type 
       filled with random numbers in the range [0,0xffffff] with bits 0xff000000 for alpha mask.  
    """
    a = 0xffffff*np.random.random(shape) + add
    return np.require(a, dtype)

#-----------------------------

def aranged_array(shape=(40,60), dtype=np.uint32) :
    """Returns numpy array of requested shape and type filling with ascending integer numbers.
    """
    arr = np.arange(size_from_shape(shape), dtype=dtype)
    arr.shape = shape
    return arr

#-----------------------------

def print_ndarr(nda, name='', first=0, last=5) :
    """Prints array attributes, title, and a few elements in a single line. 
    """    
    if nda is None : print '%s: %s' % (name, nda)
    elif isinstance(nda, tuple) : print_ndarr(np.array(nda), 'ndarray from tuple: %s' % name)
    elif isinstance(nda, list)  : print_ndarr(np.array(nda), 'ndarray from list: %s' % name)
    elif not isinstance(nda, np.ndarray) : print '%s: %s' % (name, type(nda))
    else: print '%s:  shape:%s  size:%d  dtype:%s %s...' % \
         (name, str(nda.shape), nda.size, nda.dtype, nda.flatten()[first:last])

#-----------------------------
#-----------------------------

if __name__ == '__main__':

    print_ndarr(random_exponential(), 'random_exponential')
    print_ndarr(random_standard(), 'random_standard')
    print_ndarr(random_1(), 'random_1', last=10)
    print_ndarr(random_256(), 'random_256', last=10)
    print_ndarr(random_xffffffff(), 'random_xffffffff')
    print_ndarr(random_standard(), 'random_standard')
    print_ndarr(aranged_array(), 'aranged_array')
    #print_ndarr(, '')
    print 'Test is completed'

#-----------------------------
