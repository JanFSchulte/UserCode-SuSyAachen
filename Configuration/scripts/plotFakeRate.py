#!/usr/bin/env python2.6

def getGraph(psets, var, others):
  from ROOT import TGraphAsymmErrors
  from array import array
  x, exl, exh = psets["center"].getX(var)
  y = psets["center"].getY(var, others)
  eyl = psets["lower"].getY(var, others, relativeTo = y)
  eyh = psets["upper"].getY(var, others, relativeTo = y)
  result = TGraphAsymmErrors(len(x), array("f",x), array("f",y), array("f",exl), array("f",exh), array("f",eyl), array("f",eyh))
  #result.SetName("FakeRate:%s %s:%s"%(psets["name"],var,others["name"]))
  return result

def getHashable(var, other):
  return "%s-%s"%(var,"_".join([ "%s:%s"%(k,other[k]) for k in other]))

def getRatio(num, denom, var, others):
  from ROOT import TGraphAsymmErrors, TF1
  from array import array
  from math import sqrt
  
  nX, nExl, nExh = num["center"].getX(var)
  nY = num["center"].getY(var, others)
  nEyl = num["lower"].getY(var, others, relativeTo = nY)
  nEyh = num["upper"].getY(var, others, relativeTo = nY)

  dX, dExl, dExh = denom["center"].getX(var)
  dY = denom["center"].getY(var, others)
  dEyl = denom["lower"].getY(var, others, relativeTo = dY)
  dEyh = denom["upper"].getY(var, others, relativeTo = dY)

  y = [1. - (n * 1./d) for n,d in zip(nY, dY)]
  eyl = [sqrt((1./d*eN)**2 + (n*1./d**2*eD)**2) for (n, eN, d, eD) in zip(nY, nEyl, dY, dEyl)]
  eyh = [sqrt((1./d*eN)**2 + (n*1./d**2*eD)**2) for (n, eN, d, eD) in zip(nY, nEyh, dY, dEyh)]

  correctionFit = TF1("correction %s %s"%(var, others), "[0]", 20., max(nX))
  correctionFit.SetParameter(0,y[0])
  
  graph = TGraphAsymmErrors(len(nX), array("f",nX), array("f",y), array("f",nExl), array("f",nExh), array("f",eyl), array("f",eyh))
  graph.Fit(correctionFit, "Q0R")
  
  return graph, correctionFit.GetParameter(0), correctionFit.GetParError(0)


def applyGraphOptions(cfg, graph, var, name, binName, i):
  def getAttribute(sec):
    attributeList = cfg.get(sec, name).split()
    return attributeList[i%len(attributeList)]
  graph.SetName("%s_%s"%(name, binName))
  graph.SetTitle("#splitline{%s}{%s}"%(name, binName))
  graph.SetLineColor(int(getAttribute("color")))
  graph.SetLineWidth(int(getAttribute("lineWidth")))
  graph.SetMarkerColor(int(getAttribute("color")))
  graph.SetMarkerStyle(int(getAttribute("marker")))
  graph.SetMarkerSize(float(getAttribute("markerSize")))
  graph.SetFillColor(int(getAttribute("fillColor")))
  graph.SetFillStyle(int(getAttribute("fillStyle")))

  
def main(argv=None):
  import sys, os
  from optparse import OptionParser
  from helpers import BetterConfigParser, parsePSet
  from ROOT import TMultiGraph, TCanvas, TLegend
  from Styles import tdrStyle
  tdrStyle()
  if argv == None:
        argv = sys.argv[1:]
  parser = OptionParser()
  parser.add_option("-c", "--config", dest="config", action="append", default=[],
                    help="Main configuration file. Can be given multiple times in case of split configurations. Default is default.ini")
  parser.add_option("-p", "--pSets", dest="pSets", action="append", default=[],
                    help="PSets to parse")

  parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
                    help="Verbose mode.")
  parser.add_option("-g", "--graph", action="store_true", dest="graph", default=False,
                    help="drawGraph")
  parser.add_option("-r", "--ratio", action="store_true", dest="ratios", default=False,
                    help="ration between first and all other pSets")

  parser.add_option("-s", "--show", action="store_true", dest="show", default=False,
                    help="wait after each plot has been drawn")

  (opts, args) = parser.parse_args(argv)

  if opts.config == []:
      opts.config = [ "default.ini" ]

  config = BetterConfigParser()
  config.read(opts.config)
  
  vars = config.get("general","variables").split()
  graphs = {}
  ratios = {}
  for var in vars:
    names = set([])
    graphs[var] = TMultiGraph(var, "%s;%s;%s"%(var,config.get("xTitle",var),config.get("yTitle",var)))
    ratios[var] = TMultiGraph("ratio_%s"%var, "%s;%s;%s"%(var,config.get("xTitle",var),config.get("yRatioTitle",var)))
    i = 0
    ratioTo = {}
    for pSetPath in opts.pSets:
      pSets = parsePSet(pSetPath)
      others = pSets["center"].getOthers(var)
      for other in others:
        skip = False
        for absvar in config.get("general","absvars").split():
          if absvar in other["values"] and other["values"][absvar] < 0:
             skip = True
        if skip:
          continue
        
        if opts.graph:
          graph = getGraph(pSets, var, other["values"])
          name = os.path.split(pSetPath)[1].replace("_cff.py","")
          applyGraphOptions(config, graph, var, name, other["title"], i)
          names.add(name)
          graphs[var].Add(graph)
        if opts.ratios:
          if pSetPath == opts.pSets[0]:
            ratioTo[getHashable(var, other["values"])] = pSets
          elif getHashable(var, other["values"]) in ratioTo:
            ratio, corr, corrErr = getRatio(pSets, ratioTo[getHashable(var, other["values"])],var, other["values"] )
            name = os.path.split(pSetPath)[1].replace("_cff.py","")
            print name, "$%s \pm %s$"%(corr, corrErr)
            names.add(name)
            applyGraphOptions(config, ratio, var, name, other["title"], i)
            ratios[var].Add(ratio)
          else:
            raise StandardError, "missing '%s' with '%s' in denominators! Got: %s"%(var, other["values"], ratioTo)
        i+=1
  
    if opts.graph or opts.ratios:
      c = TCanvas("canv_%s"%var,"Fakerates for %s"%var, 
        int(config.get("general","canvasSize").split("x")[0]), 
        int(config.get("general","canvasSize").split("x")[1]))
      c.Clear()
      if opts.ratios:
        ratios[var].Draw("ap")
        assert not opts.graph, "can not draw both, ratio and graphs"
      if opts.graph:
        graphs[var].Draw("ap")
        assert not opts.ratios, "can not draw both, ratio and graphs"
      legend = c.BuildLegend()
      legend.SetX1(float(config.get("legendPosition","default").split()[0]))
      legend.SetY1(float(config.get("legendPosition","default").split()[1]))
      legend.SetX2(float(config.get("legendPosition","default").split()[2]))
      legend.SetY2(float(config.get("legendPosition","default").split()[3]))
      legend.Draw()
      c.Update()
      for figFormat in config.get("general","figFormats").split():
        addition=""
        if opts.ratios: addition = "ratio_"
        c.Print(os.path.join(config.get("general","figureDir"), "%s%s_%s.%s")%(addition,var,"-".join(names),figFormat ),figFormat)
  
    if opts.show:
      raw_input("showing "+var)
    
if (__name__ == "__main__"):
  main()
