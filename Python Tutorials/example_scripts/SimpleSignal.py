import matplotlib.pyplot as plt
import numpy as np
import MDSplusTree


class SimpleSignal:
    """A signal class to provide an organizational scheme for storing
    timeseries data retrieved from an MDSplus signal node. To collect
    the data from an MDSplus signal node, provide the shot number,
    nodepath, and tree.
    """

    def __init__(self, shot, nodepath, tree = 'mst'):
        # Initialize the fields of the SimpleSignal class to be the input of 
        # the constructor.
        self.shot       = shot
        self.nodepath   = nodepath
        self.tree       = tree
        self.data       = None
        self.units      = None
        self.time       = None
        self.time_units = None
        self.error      = None
        self.name       = None
        
        # Get name of signal, last name of the nodepath
        self.name = self.nodepath.split(':')[-1]
        
        # Tree could have more data associated with a Signal, attempt to get this information,
        # but deal with the possibility that this data may not exist for a given signal.        
        
        try:
            self.data = MDSplusTree.getSignal(self.tree, self.shot, self.nodepath)
        except:
            print "{0} is not available for shot {1}".format(self.name, self.shot)
            return

        # Get the time stamps. We'll get into trouble with some of the
        # plotting methods in the Signal class if we don't have a time
        # base, so let's make one up if we can't find it.
        try: 
            self.time = MDSplusTree.getTime(self.tree, self.shot)
        except:
            self.time = np.arange(len(self.data))
 
        # Let's see if there's any more meta data like units and error bars.
        try:
            self.units = self.MDSplus.getUnits(self.tree, self.shot)
        # If there isn't antying stored in the units part of the node, we'll 
        # assume that the data are stored in terms of measured voltage.
            if not self.units:
                self.units = 'V'
        except:
            self.units = 'V'
            
        try:
            self.time_units = MDSplusTree.getTimeUnits(self.tree, self.shot)
        # If there isn't anything stored in the units part of the dim_of() node,
        # we'll assume that the data is stored in terms of measured voltage.
            if not self.time_units:
                self.time_units = 'ms'
        except:
            self.time_units = 'ms'

        # Next see whether there is any uncertainty information.
        try: 
            self.error = MDSplusTree.getError()
        except: pass
            

    def plot(self, title = None, color = None, hold = False):
        """Plot the signal vs. time using as much of the stored information
         as is available to annotate the figure."""

        if not self.time or not self.data:
            print "time and/or data is null, enter values where needed before plotting"
            return

        # Create a dictionary to pass keywords to the matplotlib plot command.
        keywords = {}
  
        # Plot data
        if color in ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']:
            keywords['color'] = color
        else:
            keywords['color'] = 'b'

        # Define the label based on the signal name including units if
        # available
        if self.name:
            keywords['label'] = self.name

        # We'll issue to plot command    
        plt.plot(self.time, self.data, **keywords)
        
        # Add a title if one is provided.
        if title:
            plt.title(title)

        # Add what annotation we can
        try: plt.xlabel(self['time_units'])
        except: pass

        try: plt.ylabel(self['units'])
        except: pass

        plt.legend()

        # Display graph
        plt.show()
