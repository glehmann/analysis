#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itk, sys, traceback
from ExtraData import ExtraData
itk.auto_not_in_place()

nucleusReader = itk.ImageFileReader.IUC3.New()
centromeresReader = itk.ImageFileReader.IUC3.New()
centromeresIntReader = itk.bioformats()

nucleusLM = itk.LabelImageToShapeLabelMapFilter.IUC3LM3.New(nucleusReader)

singleNucleusLM = itk.LabelSelectionLabelMapFilter.LM3.New(nucleusLM, InPlace=False)
binaryNucleus = itk.LabelMapToBinaryImageFilter.LM3IUC3.New(singleNucleusLM, ForegroundValue=0, BackgroundValue=255)
maurerSingleNucleus = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(binaryNucleus, UseImageSpacing=True, SquaredDistance=False)
evfSingleNucleus = itk.ErodedVolumeFractionMapImageFilter.IF3IF3.New(maurerSingleNucleus)

maskedCentromeres = itk.LabelMapMaskImageFilter.LM3IUC3.New(nucleusLM, centromeresReader)
centromeresLM = itk.LabelImageToStatisticsLabelMapFilter.IUC3IUC3LM3.New(maskedCentromeres, centromeresIntReader)

def interpolate(image, pos):
  li = itk.LinearInterpolateImageFunction.IF3D.New(image)
  return li.Evaluate( pos )

extra = ExtraData(sys.argv[1], sys.argv[2:],
          {"nucleus-ext": "-nuclei.nrrd",
           "spots-ext": "-CENP.nrrd",
           "spots-are-labeled": "true",
           "rawspots-ext": "",
           "rawspots-channel": "0"})
cf = extra.printer(file(sys.argv[1].replace("::","-")+".txt", "w"))
cf.printHeader("image", "label", "clabel", "mdist", "evf")

for (nucleusName, centromeresName, centromeresIntName) in extra.iterate(["nucleus-ext", "spots-ext", "rawspots-ext"]):
  try:
    nucleusReader(FileName=nucleusName)
    centromeresReader(FileName=centromeresName)
    centromeresIntReader(FileName=centromeresIntName, Channel=extra.getint("rawspots-channel"))
    
    for nucleus in nucleusLM()[0]:
      singleNucleusLM.SetLabel(nucleus.GetLabel())
      maskedCentromeres.SetLabel(nucleus.GetLabel())
  
      # les mesures sur chaque centromere
      evfSingleNucleus()
      for cent in centromeresLM()[0]:
        cog = cent.GetCenterOfGravity()
        evf = interpolate(evfSingleNucleus, cog)
        d = interpolate(maurerSingleNucleus, cog)
        cf.printData(nucleusName, nucleus.GetLabel(), cent.GetLabel(), d, evf)

      # flush every nucleus to be able to check the result quickly
      cf.flush()

  except Exception, e:
    print >> sys.stderr, "Error with", nucleusName
    print >> sys.stderr, e
    traceback.print_exc(file=sys.stderr)
    centromeresLM.ResetPipeline()
    evfSingleNucleus.ResetPipeline()

