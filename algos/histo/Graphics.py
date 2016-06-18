#------------------------------
import numpy as np

import matplotlib
#if matplotlib.get_backend() != 'Qt4Agg' : matplotlib.use('Qt4Agg')

import matplotlib.pyplot  as plt
#import matplotlib.lines   as lines
#import matplotlib.patches as patches

from CalibManager.PlotImgSpeWidget import add_stat_text

#------------------------------
#class Storage :
#    def __init__(self) :
#        pass
#
#------------------------------
#store = Storage() # singleton
#------------------------------

#------------------------------

def figure(figsize=(13,12), title='Image', dpi=80, facecolor='w', edgecolor='w', frameon=True,
           move=None\
          ) :
    """ Creates and returns figure
    """
    fig = plt.figure(figsize=figsize,\
                     dpi=dpi,\
                     facecolor=facecolor,\
                     edgecolor=edgecolor,\
                     frameon=frameon)
    fig.canvas.set_window_title(title)
    if move is not None : move_fig(fig, x0=move[0], y0=move[1])
    return fig

#------------------------------

def move_fig(fig, x0=200, y0=100) :
    fig.canvas.manager.window.geometry('+%d+%d' % (x0, y0))

#------------------------------

def move(x0=200,y0=100) :
    plt.get_current_fig_manager().window.geometry('+%d+%d' % (x0, y0))

#------------------------------

def add_axes(fig, axwin=(0.05, 0.03, 0.87, 0.93)) :
    """Add axes to figure from input list of windows.
    """
    axes = fig.add_axes(axwin)
    return axes

#------------------------------

def set_win_title(fig, titwin='Image') :
    fig.canvas.set_window_title(titwin)

#------------------------------

def add_title_labels_to_axes(axes, title=None, xlabel=None, ylabel=None, fslab=14, fstit=20, color='k') :
    if title  is not None : axes.set_title(title, color=color, fontsize=fstit)
    if xlabel is not None : axes.set_xlabel(xlabel, fontsize=fslab)
    if ylabel is not None : axes.set_ylabel(ylabel, fontsize=fslab)

#------------------------------

def show(mode=None) :
    if mode is None : plt.ioff() # hold contraol at show() (connect to keyboard for controllable re-drawing)
    else            : plt.ion()  # do not hold control
    plt.show()

#------------------------------

def draw(fig) :
    plt.draw()

#------------------------------

def draw_fig(fig) :
    fig.canvas.draw()

#------------------------------

def save_plt(fname='img.png', verb=True) :
    if verb : print 'Save plot in file: %s' % fname 
    plt.savefig(fname)

#------------------------------

def save_fig(fig, fname='img.png', verb=True) :
    if verb : print 'Save figure in file: %s' % fname 
    fig.savefig(fname)

#------------------------------

def hist(axhi, arr, bins=None, amp_range=None, weights=None, color=None, log=False) :
    """Makes historgam from input array of values (arr), which are sorted in number of bins (bins) in the range (amp_range=(amin,amax))
    """
    #axhi.cla()
    hi = axhi.hist(arr.flatten(), bins=bins, range=amp_range, weights=weights, color=color, log=log) #, log=logYIsOn)
    if amp_range is not None : axhi.set_xlim(amp_range) # axhi.set_autoscale_on(False) # suppress autoscailing
    wei, bins, patches = hi
    add_stat_text(axhi, wei, bins)
    return hi

#------------------------------

def imshow(axim, img, amp_range=None, extent=None,\
           interpolation='nearest', aspect='auto', origin='upper',\
           orientation='horizontal', cmap='jet') :
    """
    extent - list of four image physical limits for labeling,
    cmap: 'gray_r'
    #axim.cla()
    """
    imsh = axim.imshow(img, interpolation=interpolation, aspect=aspect, origin=origin, extent=extent, cmap=cmap)
    if amp_range is not None : imsh.set_clim(amp_range[0],amp_range[1])
    return imsh

#------------------------------

def colorbar(fig, imsh, axcb, orientation='vertical', amp_range=None) :
    """
    orientation = 'horizontal'
    amp_range = (-10,50)
    """
    if amp_range is not None : imsh.set_clim(amp_range[0],amp_range[1])
    cbar = fig.colorbar(imsh, cax=axcb, orientation=orientation)
    return cbar

#------------------------------

def imshow_cbar(fig, axim, axcb, img, amin=None, amax=None, extent=None,\
                interpolation='nearest', aspect='auto', origin='upper',\
                orientation='horizontal', cmap='jet') :
    """
    extent - list of four image physical limits for labeling,
    cmap: 'gray_r'
    #axim.cla()
    """
    axim.cla()
    imsh = axim.imshow(img, interpolation=interpolation, aspect=aspect, origin=origin, extent=extent, cmap=cmap)
    cbar = fig.colorbar(imsh, cax=axcb, orientation=orientation)
    return imsh, cbar

#------------------------------
#------------------------------
#------------------------------
#------------------------------

def test01() :
    """ imshow
    """
    from pyimgalgos.NDArrGenerators import random_standard

    img = random_standard(shape=(40,60), mu=200, sigma=25)
    fig = figure(figsize=(6,5), title='Test imshow', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))    
    axim = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))
    imsh = imshow(axim, img, amp_range=None, extent=None,\
           interpolation='nearest', aspect='auto', origin='upper',\
           orientation='horizontal', cmap='jet') 
    show()

#------------------------------

def test02() :
    """ hist
    """
    from pyimgalgos.NDArrGenerators import random_standard

    mu, sigma = 200, 25
    arr = random_standard((500,), mu, sigma)
    fig = figure(figsize=(6,5), title='Test hist', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))    
    axhi = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))
    his = hist(axhi, arr, bins=100, amp_range=(mu-6*sigma,mu+6*sigma), weights=None, color=None, log=False)
    show()

#------------------------------

def test03() :
    """ Update image in the event loop
    """
    from pyimgalgos.NDArrGenerators import random_standard

    mu, sigma = 200, 25
    fig = figure(figsize=(6,5), title='Test hist', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))
    axim = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))

    imsh = None

    for i in range(100) :
       img = random_standard((1000,1000), mu, sigma)
       #axim.cla()
       set_win_title(fig, 'Event %d' % i)

       if imsh is None :
           imsh = imshow(axim, img, amp_range=None, extent=None,\
                  interpolation='nearest', aspect='auto', origin='upper',\
                  orientation='horizontal', cmap='jet') 
       else :
           imsh.set_data(img)

       show(mode=1)  # !!!!!!!!!!       
       draw_fig(fig) # !!!!!!!!!!
    show()

#------------------------------

def test04() :
    """ Update histogram in the event loop
    """
    from pyimgalgos.NDArrGenerators import random_standard

    mu, sigma = 200, 25
    fig = figure(figsize=(6,5), title='Test hist', dpi=80, facecolor='w', edgecolor='w', frameon=True, move=(100,10))
    axhi = add_axes(fig, axwin=(0.10, 0.08, 0.85, 0.88))

    for i in range(10) :
       arr = random_standard((500,), mu, sigma)
       axhi.cla()
       set_win_title(fig, 'Event %d' % i)
       his = hist(axhi, arr, bins=100, amp_range=(mu-6*sigma,mu+6*sigma), weights=None, color=None, log=False)

       draw(fig)    # !!!!!!!!!!
       show(mode=1) # !!!!!!!!!!
    show()

#------------------------------
#------------------------------
#------------------------------
#------------------------------

def test_selected() :

    import sys

    if len(sys.argv)==1   :
        print 'Use command > python %s <test-number [1-5]>' % sys.argv[0]
        sys.exit ('Add <test-number> in command line...')

    elif sys.argv[1]=='1' : test01()
    elif sys.argv[1]=='2' : test02()
    elif sys.argv[1]=='3' : test03()
    elif sys.argv[1]=='4' : test04()
    else :
        print 'Non-expected arguments: sys.argv=', sys.argv
        sys.exit ('Check input parameters')

#------------------------------

def test_all() :
    test01()
    test02()

#------------------------------

if __name__ == "__main__" :

    test_selected()
    #test_all()
    print 'End of test'

#------------------------------

