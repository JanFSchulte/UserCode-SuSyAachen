#!/usr/bin/env python

def parsePSet(path):
  import cmsDummies
  result = {}
  file = open(path, "r")
  rawPSets = file.read()
  file.close()
  execGlobals = {"cms":cmsDummies}
  exec(rawPSets, execGlobals )
  for psetName in filter(lambda x: x not in ["cms", "__builtins__"], execGlobals.keys()):
    if "Center" in psetName:
      result["name"] = psetName.split("Center")[0]
      result["center"] = execGlobals[psetName]
    if "Lower" in psetName:
      result["lower"] = execGlobals[psetName]
    if "Upper" in psetName:
      result["upper"] = execGlobals[psetName]
  return result

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
  from helpers import BetterConfigParser
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
  parser.add_option("-s", "--show", action="store_true", dest="show", default=False,
                    help="wait after each plot has been drawn")

  (opts, args) = parser.parse_args(argv)

  if opts.config == []:
      opts.config = [ "default.ini" ]

  config = BetterConfigParser()
  config.read(opts.config)
  
  vars = config.get("general","variables").split()
  graphs = {}
  for var in vars:
    names = set([])
    graphs[var] = TMultiGraph(var, "%s;%s;%s"%(var,config.get("xTitle",var),config.get("yTitle",var)))
    i = 0
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
          i+=1
  
    if opts.graph:
      c = TCanvas("canv_%s"%var,"Fakerates for %s"%var, 
        int(config.get("general","canvasSize").split("x")[0]), 
        int(config.get("general","canvasSize").split("x")[1]))
      graphs[var].Draw("ap")
      legend = c.BuildLegend()
      legend.SetX1(float(config.get("legendPosition","default").split()[0]))
      legend.SetY1(float(config.get("legendPosition","default").split()[1]))
      legend.SetX2(float(config.get("legendPosition","default").split()[2]))
      legend.SetY2(float(config.get("legendPosition","default").split()[3]))
      legend.Draw()
      c.Update()
      for figFormat in config.get("general","figFormats").split():
        c.Print(os.path.join(config.get("general","figureDir"), "%s_%s.%s")%(var,"-".join(names),figFormat ),figFormat)
  
    if opts.show:
      raw_input("showing "+var)
    
if (__name__ == "__main__"):
  main()