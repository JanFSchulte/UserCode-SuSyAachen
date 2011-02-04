class VPSet:
  def __init__(self, *args):
    self.binBoundaries = {}
    self.bins = args
    for bin in self.bins:
      for opt in dir(bin):
        if opt.endswith("Min") or opt.endswith("Max"):
          if not opt[:-3] in self.binBoundaries:
            self.binBoundaries[opt[:-3]] = set([])
          self.binBoundaries[opt[:-3]].add( getattr( bin, opt) )
    for bound in self.binBoundaries:
      self.binBoundaries[bound] = sorted(self.binBoundaries[bound])
      
  def getX(self, var):
    from math import fabs
    x = []
    low = []
    high = []
    for i in range(len(self.binBoundaries[var])-1):
      x.append(0.5 * (self.binBoundaries[var][i+1] + self.binBoundaries[var][i]))
      low.append(fabs( self.binBoundaries[var][i] - x[-1]))
      high.append(fabs(self.binBoundaries[var][i+1]-x[-1]))
    return (x, low, high)
      
  def getY(self, var, others, relativeTo = None):
    from math import fabs
    result = []
    x, exl, exh = self.getX(var)
    for i in range(len(x)):
      frPoint = {var: x[i]}
      frPoint.update(others)
      y = self.fakeRate(frPoint)
      if not relativeTo == None:
        y = fabs(y - relativeTo[i])
        
      result.append(y)
    return result
  
  def fakeRate(self, point):
    result = None
    for bin in self.bins:
      useBin = True
      for var, val in point.iteritems():
 #       print not "%sMin"%var in dir(bin), getattr(bin, "%sMin"%var), "<=",val
 #       print not "%sMax"%var in dir(bin), getattr(bin, "%sMax"%var),">", val
        useBin &= (not "%sMin"%var in dir(bin)) or getattr(bin, "%sMin"%var) <= val
        useBin &= (not "%sMax"%var in dir(bin)) or getattr(bin, "%sMax"%var) > val
      if useBin: result = bin.weight
#      print "--",useBin,"--"
    assert result != None, "failed to calculate fakerate for '%s' in '%s'"%(point, [i.__dict__ for i in self.bins])
    return result
    
  def getOthers(self, chosenVar):
    result = [{"values":{}, "name":""}]
    for otherVar in filter(lambda x: not x == chosenVar, self.binBoundaries.keys()):
      x, exl, exh = self.getX(otherVar)
      tmp = []
      for previous in result:
        for i in range(len(x)):
          values = {otherVar: x[i]}
          values.update(previous["values"])
          tmp.append({ "values": values,"name":"%s #epsilon [%s, %s) %s"%(otherVar, x[i]-exl[i], x[i]+exh[i], previous["name"])})
      result = tmp
    return result
    
class PSet:
  def __init__(self, **kwargs):
    #print self.__name__
    for arg in kwargs:
      self.__dict__[arg] = kwargs[arg].value
    
class double:
  def __init__(self, value):
    #print self.__name__
     self.value = value
  

