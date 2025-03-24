import numpy as np
from . import Source
from .. import util
import time

class segmentedDike(Source):
    """
    A class used to represent a segmented dike model composed of multiple connected patches.
    
    attributes : 
    nl : number of segments in the lenghth direction
    nw : number of segments in the width direction

    parameters : 
    xcen : x-coordinate for the center of the dike
    ycen : y-coordinate for the center of the dike
    depth : depth of the center of the dike
    length : length of the dike (horizontal derection)
    width : width of the dike (vertical direction)
    strike : orientation of the dike clockwise from north
    openning: opennings of segments
    """

    def __init__(self, data, nl=None, nw=None):
        self.nl     = 1 if nl     is None else nl
        self.nw     = 1 if nw     is None else nw

        super().__init__(data)
    
    def get_source_id(self):
        return "segmentedDike"
    
    def bayesian_steps(self):
        steps = 1010000  
        burnin = 10000   
        thin = 1000      
        return steps, burnin, thin
        
    def time_dependent(self):
        return False
    
    def set_parnames(self):
        """
        Function defining the names for the parameters in the model.
        """
        self.parameters = ("xcen", "ycen", "depth", "length", "width", "strike")
        for i in range(self.nl*self.nw):
            self.parameters += ('open'+str(i),)    

    def get_patches(self, xcen, ycen, depth, length, width, strike):
        """
        list of tuples (x, y, d), center coordinates of patches
        return order ：vertically from up to down, and within the same row from left to right
        """
        strike_rad = np.deg2rad(strike)
    
        patch_length = length / self.nl
        patch_width  = width  / self.nw
    
        patches = []
        # row, width direction
        for j in range(nw):
            patch_depth = depth - width/2 + (j + 0.5) * patch_width
            # row/length direction
            for i in range(nl):
                u = -length/2 + (i + 0.5) * patch_length
                offset_x = u * np.sin(strike_rad)
                offset_y = u * np.cos(strike_rad)
                patch_x = xcen + offset_x
                patch_y = ycen + offset_y
                patches.append((patch_x, patch_y, patch_depth))
        return patches

    def get_greens(self, xcen, ycen, depth, length, width, strike):
        x = self.data.xs
        y = self.data.ys

        dat1 = Gnss()
        dat1.add_xs(x)
        dat1.add_ys(y)
        dat1.add_data(x*0,x*0,x*0)

        oki = Okada(dat1)
        oki.set_type('open')

        patches = self.get_patches()
        xo = [xcen,ycen,depth,length,width,1,strike,90]
        defo=oki.forward(xo)
        G = np.zeros( (len(defo),nl*nw) )
        patch_length = length / self.nl
        patch_width  = width  / self.nw
        for i, patch in enumerate(patches):
            xp = [ patch[0],patch[1],patch[2],patch_length,patch_width,1,strike,90]
            defo = oki.forward(xp)
            G[:,i] = defo
        return G
    
    def model(self, x, y, xcen, ycen, depth, length, width, strike, *ops):
        G = self.get_greens(xcen,ycen,depth,length,width,strike)
        opennings = np.array(ops)
        data = G@opennings

        ux=data[0:len(x)]
        uy=data[len(x):2*len(x)]
        uz=data[2*len(x):3*len(x)]
        return ux,uy,uz
    
