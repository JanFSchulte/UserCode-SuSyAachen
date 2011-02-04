'''
Created on 18.08.2009

@author: heron
'''
#-------------  need to put this in its own package some time
import ConfigParser
import re
class BetterConfigParser(ConfigParser.ConfigParser):

    def __init__(self, defaults=None):#in python 2.6, dict_type=dict):
        self.__runtimeRepMap = {}
        ConfigParser.ConfigParser.__init__(self, defaults)#, dict_type)            

    def get(self, section, option, raw = True, default = None, sectionwideReplace = True):
        result = None
        matchedName = None
        for name, item in self.items(section, raw=raw):
            if re.match(name+"$", option) != None:
              assert matchedName == None, "option '%s' matches '%s' and '%s'"%(option, matchedName, name)
              matchedName = name
              result = item
        if matchedName == None:
            if not default is None:
                result = default
            if self.has_option(section, "default"):
                result = self.get(section, "default", raw=raw)
#        result = ConfigParser.ConfigParser.get(self, section, matchedName, raw)
        if sectionwideReplace:
            result = "\n".join([self.__replaceSectionwideTemplates(l) for l in result.splitlines()])
        #not sure if I have created a chicken egg problem here
        for entry in self.__runtimeRepMap:
            result = result.replace(".oO[%s]Oo." % entry, self.__runtimeRepMap[entry])
        return result
    
# need to look this syntax up once I am back online...
#    def items(self, **opts):
#        opts["raw"] = True
#        return ConfigParser.ConfigParser.itmes(self, **opts)
    
    def optionxform(self, optionstr):
        '''enable case sensitive options in .ini files'''
        return optionstr

    def write(self, fileObject):
        for section in self.sections():
            for option in self.options(section):
                self.set(section, option, self.get(section,option, raw=True, sectionwideReplace=False))
        ConfigParser.ConfigParser.write(self, fileObject)

    def addRuntimeReplacement(self, entry, replacement):
        if entry in self.__runtimeRepMap:
            raise StandardError, "'%s' is already in the runtime replacements. You are not allowed to modify!" % entry
        self.__runtimeRepMap[entry] = replacement

    def __replaceSectionwideTemplates(self, data):
        '''replace <section|option> with get(section,option) recursivly'''
        result = data
        findExpression = re.compile("((.*)\<([^\>]*)\|([^\>]*)\>(.*))*")
        groups = findExpression.search(data).groups()
        if not self.has_section(groups[2]) and not groups[2] == None:
            print "WARNING: skipping '%s' in '%s' because section could not be found"%(groups[2],data)
        if not groups == (None, None, None, None, None) and self.has_section(groups[2]): # expression not matched
            result = self.__replaceSectionwideTemplates(groups[1])
            
            matchedName = None
            for name, item in self.items(groups[2], raw=True):
                if re.match(name+"$", groups[3]) != None:
                    if item == "passedString": result +=  self.__replaceSectionwideTemplates(groups[3])
                    else: result +=  self.__replaceSectionwideTemplates(item)
                    assert matchedName == None, " '<%s|%s>' did not only match '%s' but also '%s'"%(groups[2],groups[3], name, matchedName)
                    matchedName = name
                    
            if matchedName == None:
                if self.has_option(groups[2],"default"):
                    if self.get(groups[2], "default") == "passedString":
                        result +=  self.__replaceSectionwideTemplates( groups[3] )
                    else:
                        result +=  self.__replaceSectionwideTemplates( self.get(groups[2], "default") )
                else:
                    #rethrow exception for the poor
                    result += self.get(groups[2], groups[3])
            result += self.__replaceSectionwideTemplates(groups[4])
        return result

class Section:
    def __init__(self, config, name, defaults = {}):
        self._config = config
        self._name = name
        self._defaults = defaults
    
    def setDefaults(self, defaults):
        self._defaults.update(defaults)
    
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[ name ]
        else:
            assert self._name in self._config.sections(), "Could not find section '%s' while getting '%s'"%(self._name, name)
            self._defaults.update( eval(self._config.get(self._name, "defaults", default = "{}") ) )
            default = None
            if name in self._defaults: default = self._defaults[name]
            return self._config.get(self._name, name, default = default)
    
