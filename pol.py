#!/usr/bin/env python
# -*- coding: utf-8 -*-

import itk, random, sys, re, traceback
from ExtraData import ExtraData
itk.auto_not_in_place()
random.seed()

def computeRegion(lo, img):
   img = itk.output(img)
   import math
   transform = lo.GetBinaryPrincipalAxesToPhysicalAxesTransform()
   invTransform = itk.AffineTransform.D3.New()
   transform.GetInverse( invTransform )
   region = lo.GetRegion()
   outputSpacing = [min(itk.spacing(img))]*3
   dummyImg = itk.Image.UC3.New()
   dummyImg.SetSpacing( outputSpacing )
   idx = region.GetIndex()
   size = region.GetSize()
   minPoint = [4294967295L]*3
   maxPoint = [-4294967295L]*3
   for x in [idx[0], idx[0]+size[0]]:
       for y in [idx[1], idx[1]+size[1]]:
           for z in [idx[2], idx[2]+size[2]]:
               inputPoint = img.TransformIndexToPhysicalPoint( [x,y,z] )
               outputPoint = invTransform.TransformPoint( inputPoint )
               for i in range(0,3):
                   if minPoint[i] > outputPoint[i]:
                       minPoint[i] = outputPoint[i]
                   if maxPoint[i] < outputPoint[i]:
                       maxPoint[i] = outputPoint[i]
   minIdx = dummyImg.TransformPhysicalPointToIndex(minPoint)
   maxIdx = dummyImg.TransformPhysicalPointToIndex(maxPoint)
   outputIdx = []
   outputSize = []
   for i in range(0,3):
       outputIdx.append(int(minIdx[i]))
       outputSize.append( int(maxIdx[i] - minIdx[i]) )
   uabc3 = 1.0
   for v in itk.spacing(img):
     uabc3 *= v
   for v in lo.GetEquivalentEllipsoidSize():
     uabc3 *= v
   uabc3 = math.pow(uabc3, 1/3.0)
   spacing2 = [uabc3 / v for v in lo.GetEquivalentEllipsoidSize()]
   return (outputIdx, outputSize, outputSpacing, spacing2)

nucleusReader = itk.ImageFileReader.IUC3.New()
centromeresReader = itk.ImageFileReader.IUC3.New()
centromeresIntReader = itk.ImageFileReader.IUC3.New()

# unflatten the nucleus and the images associated to it
li2lm = itk.LabelImageToShapeLabelMapFilter.IUC3LM3.New(nucleusReader)
nni = itk.NearestNeighborInterpolateImageFunction.IUC3D.New()
resampleNucleus = itk.ResampleImageFilter.IUC3IUC3.New(nucleusReader, Interpolator=nni)
resampleCentromeres = itk.ResampleImageFilter.IUC3IUC3.New(centromeresReader, Interpolator=nni)
resampleCentromeresInt = itk.ResampleImageFilter.IUC3IUC3.New(centromeresIntReader)
unflattenNucleus = itk.ChangeInformationImageFilter.IUC3.New(resampleNucleus, ChangeSpacing=True)
unflattenCentromeres = itk.ChangeInformationImageFilter.IUC3.New(resampleCentromeres, ChangeSpacing=True)
unflattenCentromeresInt = itk.ChangeInformationImageFilter.IUC3.New(resampleCentromeresInt, ChangeSpacing=True)

nucleusLM = itk.LabelImageToShapeLabelMapFilter.IUC3LM3.New(unflattenNucleus)

maskedCentromeres = itk.LabelMapMaskImageFilter.LM3IUC3.New(nucleusLM, unflattenCentromeres)
centromeresLM = itk.LabelImageToStatisticsLabelMapFilter.IUC3IUC3LM3.New(maskedCentromeres, unflattenCentromeresInt)
posCentromeresLM = itk.StatisticsPositionLabelMapFilter.LM3.New(centromeresLM, Attribute="CenterOfGravity")
aggrCentromeresLM = itk.AggregateLabelMapFilter.LM3.New(posCentromeresLM, InPlace=False)
shapeAggrCentromeresLM = itk.ShapeLabelMapFilter.LM3.New(aggrCentromeresLM)

simCentromeresLM = itk.LabelMap._3.New()
shapeSimCentromeresLM = itk.ShapeLabelMapFilter.LM3.New(simCentromeresLM)

simCentromeresLM2 = itk.LabelMap._3.New()
shapeSimCentromeresLM2 = itk.ShapeLabelMapFilter.LM3.New(simCentromeresLM2)

extra = ExtraData(defaults={"nucleus-ext": "-nuclei.nrrd",
           "spots-ext": "-CENP.nrrd",
           "spots-are-labeled": "true",
           "rawspots-ext": "",
           "rawspots-channel": "0",
           "number-of-simulations": "500"})

nf = extra.printer(file(sys.argv[1].replace("::","-")+".txt", "w"))
nf.printHeader("image", "label", "size", "centromeres", "flatness", "cdist", "rdist", "ncdist", "nrdist", "pvalue")

for (nucleusName, centromeresName, centromeresIntName) in extra.iterate(["nucleus-ext", "spots-ext", "rawspots-ext"]):
  try:
    nbSims = extra.getint("number-of-simulations")
    nucleusReader(FileName=nucleusName)
    centromeresReader(FileName=centromeresName)
    centromeresIntReader(FileName=centromeresIntName)

    for flatNucleus in li2lm()[0]:    
      idx, size, spacing, spacing2 = computeRegion(flatNucleus, nucleusReader)
      transform = flatNucleus.GetBinaryPrincipalAxesToPhysicalAxesTransform()
      resampleNucleus(Transform=transform, OutputStartIndex=idx, Size=size, OutputSpacing=spacing)
      resampleCentromeres(Transform=transform, OutputStartIndex=idx, Size=size, OutputSpacing=spacing)
      resampleCentromeresInt(Transform=transform, OutputStartIndex=idx, Size=size, OutputSpacing=spacing)
      unflattenNucleus(OutputSpacing=spacing2)
      unflattenCentromeres(OutputSpacing=spacing2)
      unflattenCentromeresInt(OutputSpacing=spacing2)

      nucleus = nucleusLM()[0].GetLabelObject(flatNucleus.GetLabel())
      maskedCentromeres.SetLabel(nucleus.GetLabel())

      nucleusSize = nucleus.GetSize()
      npos = nucleus.GetCentroid()

      # une simulation pour une analyse globale dans R
      simCentromeresLM.SetRegions(itk.region(nucleusReader))
      simCentromeresLM.SetSpacing(itk.spacing(nucleusReader))
      simCentromeresLM.ClearLabels()

      for i in range(centromeresLM()[0].GetNumberOfLabelObjects()):
        idx = nucleus.GetIndex(random.randint(0,  nucleusSize-1))
        simCentromeresLM.SetPixel(idx, 1)
      simCentromeresLM.GetLabelObject(1).Optimize()
      
      cpos = shapeAggrCentromeresLM()[0].GetNthLabelObject(0).GetCentroid()
      cdist = (npos - cpos).GetNorm()

      spos = shapeSimCentromeresLM()[0].GetNthLabelObject(0).GetCentroid()
      sdist = (npos - spos).GetNorm()
      
      # realisation dun test par noyau avec plusieurs simulations
      simCentromeresLM2.SetRegions(itk.region(nucleusReader))
      simCentromeresLM2.SetSpacing(itk.spacing(nucleusReader))

      sdist2s = []
      for s in range(nbSims):
        simCentromeresLM2.ClearLabels()
        for i in range(centromeresLM()[0].GetNumberOfLabelObjects()):
          idx = nucleus.GetIndex(random.randint(0,  nucleusSize-1))
          simCentromeresLM2.SetPixel(idx, 1)

        simCentromeresLM2.GetLabelObject(1).Optimize()
        spos2 = shapeSimCentromeresLM2()[0].GetNthLabelObject(0).GetCentroid()
        sdist2 = (npos - spos2).GetNorm()
        sdist2s.append(sdist2)
      
      sdist2s.sort()
      sdist2s.reverse()
      pos=0
      ok=False
      while pos<nbSims and not ok:
        if cdist > sdist2s[pos]:
          ok=True
        else:
          pos += 1
      pvalue = float(pos)/nbSims

      # pour normaliser les distances dans la première méthode
      # msdist2 = sum(sdist2s)/nbSims
      msdist2 = nucleus.GetEquivalentRadius()
  
      nf.printData(nucleusName, nucleus.GetLabel(), nucleus.GetSize(), centromeresLM()[0].GetNumberOfLabelObjects(), nucleus.GetBinaryFlatness(), cdist, sdist, cdist/msdist2, sdist/msdist2, pvalue)
      nf.flush()

  except Exception, e:
    print >> sys.stderr, "Error with", nucleusName
    print >> sys.stderr, e
    traceback.print_exc(file=sys.stderr)
    shapeAggrCentromeresLM.ResetPipeline()

