####################################
# CLASS TO COMPUTE RUNNING STATISTICS
"""
BASED on Welford's algorithim to incrementally
compute sample variance
"""
###################################
import numpy as np

class OnlineVariance(object):

      def __init__(self,iterable=None,ddof=1):
	 self.ddof,self.n,self.mean,self.M2=ddof, 0,0.0,0.0
	 if iterable is not None:
            for datum in iterable:
		self.include(datum)


      def include(self,datum):
	  self.n+=1
	  self.delta =datum-self.mean
	  self.mean+=self.delta/self.n
	  self.M2+=self.delta*(datum-self.mean)
	  self.variance = self.M2 / (self.n - self.ddof)

      @property
      def std(self):
	  return np.sqrt(self.variance)




class OnlineSkewKurt(object):
      def __init__(self,iterable=None):
	  self.n = 0
	  self.mean =0
	  self.M2 = 0
	  self.M3=0
	  self.M4=0
	  if iterable is not None:
             for datum in iterable:
                 self.include(datum)
      def include(self,datum):
          self.n1 = self.n
	  self.n +=1
	  self.delta = datum-self.mean
	  self.delta_n = self.delta / self.n
	  self.delta_n2 = self.delta_n * self.delta_n
	  self.term1 = self.delta * self.delta_n * self.n1
	  self.mean+=self.delta_n
	  self.M4 = self.M4 + self.term1 * self.delta_n2 * (self.n*self.n - 3*self.n + 3)\
			  + 6*self.delta_n2 * self.M2 -4*self.delta_n*self.M3
	  self.M3+=self.term1*self.delta_n*(self.n-2)-3*self.delta_n*self.M2
	  self.M2+=self.term1

      def kurt_skew(self):
	  self.kurtosis = (self.n*self.M4) / (self.M2*self.M2) -3
	  self.skew    = (np.sqrt(self.n)*self.M3) / (self.M2**1.5)


