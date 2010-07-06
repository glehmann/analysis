# -*- coding: utf-8 -*-
import ConfigParser, os, glob, sys, commands, time

class ExtraData:
  attrs=None
  paths=None
  def __init__(self, section=None, paths=[], defaults={}):
    if section == None:
      if len(sys.argv) >= 2:
        section = sys.argv[1]
      elif defaults.has_key("section"):
        section = defaults["section"]
      else:
        section = "analysis"
    self.section = section
    if paths == []:
      # lets use the paths provided on the command line
      paths = sys.argv[2:]
    if paths == []:
      # we have to find the paths automatically
      paths = commands.getoutput("find . -name extradata.txt").split()
    if paths == []:
      # still no path - lets use the working dir
      paths = ["."]
    self.paths = paths
    self.defaults = defaults
    self.read(paths[0])

  def read(self, path):
    if os.path.isdir(path):
       self.dirname = path
       self.config = ConfigParser.RawConfigParser()
       self.config.add_section("experiment")
    else:
      self.dirname = os.path.dirname(path)
      if self.dirname == "":
        self.dirname = "."
      self.config = ConfigParser.RawConfigParser()
      self.config.read(path)
      if self.config.has_section("include") and self.config.has_option("include", "file"):
        self.config.read(self.dirname+os.sep+self.config.get("include", "file"))
        # reread the main file to overide the ones in the included file
        self.config.read(path)
    if self.attrs == None:
      self.attrs = tuple(self.config.options("experiment"))
    
  def iterate(self, attribNames, quiet=False):
    for p in self.paths:
      self.read(p)
      exts = [self.get(a) for a in attribNames]
      ref, others = exts[0], exts[1:]
      fs = glob.glob(self.dirname+os.sep+'*'+ref)
      if fs == [] and not quiet:
        print >> sys.stderr, "Warning: no such file:", self.dirname+os.sep+'*'+ref
      for f in fs:
        base = f[:-len(ref)]
        r = [f]
        for e in others:
          if os.path.exists(base+e):
            r.append(base+e)
          elif not quiet:
            print >> sys.stderr, "Warning: %s is missing" % (base+e)
        if len(r) == len(exts):
          if not quiet:
            print >> sys.stderr, "Processing", f
          if len(r) == 1:
            yield r[0]
          else:
            yield tuple(r)
          
  def get(self, k):
    if self.config.has_option(self.section, k):
      return self.config.get(self.section, k)
    return self.defaults[k]

  def getint(self, k):
    return int(self.get(k))

  def getfloat(self, k):
    return float(self.get(k))

  def getbool(self, k):
    return self.get(k).lower() == "true" or self.get(k).lower() == "on" or self.get(k) == "1"

  def printer(self, out=sys.stdout):
    return Printer(self, out)

  def basename(self):
    return self.section.replace("::","-")

class Printer:
  def __init__(self, extra, out=sys.stdout):
    self.extra = extra
    self.out = out
    if isinstance(self.out, str):
      if os.path.exists(self.out):
        t = time.localtime(os.path.getmtime(self.out))
        ext = "-%s-%s-%s_%sh%s.bak" % (t.tm_year, t.tm_mon, t.tm_mday, t.tm_hour, t.tm_min)
        os.rename(self.out, self.out+ext)
        os.utime(self.out+ext, None)
      self.out = file(self.out, "w")
    self.lenght = 0

  def printHeader(self, *extendedHeaders):
    headers = self.extra.attrs+extendedHeaders
    self.length = len(headers)
    print >> self.out, '"%s"' % '" "'.join(headers)
    self.out.flush()

  def printData(self, *extendedData):
    data = [self.extra.config.get("experiment", k) for k in self.extra.attrs]+[str(d) for d in extendedData]
    if len(data) != self.length:
      raise Exception("Data and header length mismatch - got %i, expected %i." % (len(data), self.length))
    print >> self.out, ' '.join(data)
    self.out.flush()
  
  def close(self):
    self.out.close()
  
  def flush(self):
    self.out.flush()
