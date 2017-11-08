import random
import numpy as np
class MDSplusTree:

	def __init__(self, units = None, time = None):
		random.seed(7)
		self.trees = ['mst', 'plasma', 'physics']
		self.ydata = dict() 
		self.xdata = dict()
		self.units = units
		self.time  = time
		
		self.xdata['mst']     = [i for i in range(-50, 50)]
		self.xdata['plasma']  = [i for i in range(100)]
		self.xdata['physics'] = [i for i in range(0,200,2)]
		
		for i in range(999):
			if i % 3 == 0:
				self.ydata[i] = self._normal
			elif i % 3 == 1:
				self.ydata[i] = self._rising
			else:
				self.ydata[i] = self._random
	
	@property
	def _normal(self):
		scale = random.random()
		first  = [scale*float(i) for i in range(0,50)]
		second = [scale*float(i) for i in range(49,-1,-1)]
		first.extend(second)
		return first
		
	@property
	def _rising(self):
		scale = random.random()
		return [scale*float(i) for i in range(100)]
		
	@property
	def _random(self):
		scale = random.random()
		return list(scale*np.random.rand(100))
		
	def getSignal(self, tree, shot):
		return [self.xdata[tree], self.ydata[shot]]
		
	def getUnits(self, tree, shot):
		return self.units
		
	def getTime(self, tree, shot):
		return self.time