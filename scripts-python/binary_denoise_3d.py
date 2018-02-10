import itk
import sys
import os
if len(sys.argv) != 4:
  print("Usage: " + sys.argv[0] + "inputImage outputFOLDER physicalsize_lambda")
  sys.exit(1)
print("binary_denoise_3d %s" % sys.argv[1])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename = os.path.splitext(basename)[0]
output_folder = sys.argv[2]
physical_size_lambda = sys.argv[3]
output_extension = "nrrd"
output_filename = os.path.join(output_folder, filename +
                               "_Size" + physical_size_lambda +
                               "." + output_extension)
physical_size_lambda = float(physical_size_lambda)

PixelType = itk.UC
Dimension = 3
ImageType = itk.Image[PixelType, Dimension]
LabelObjectType = itk.StatisticsLabelObject[itk.UL, Dimension]
LabelMapType = itk.LabelMap[LabelObjectType]
reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()

# min_max_calculator = itk.MinimumMaximumImageCalculator.New(image)
# min_max_calculator.Compute()
# max_value = min_max_calculator.GetMaximum()
# min_value = min_max_calculator.GetMinimum()
# watershed_level_normalized = (max_value - min_value)*watershed_level
# DistancePixelType = itk.F
# DistanceImageType = itk.Image[ DistancePixelType, Dimension ]
# maurer = itk.SignedMaurerDistanceMapImageFilter[ImageType, DistanceImageType].New(image)
# watershed = itk.MorphologicalWatershedImageFilter[DistanceImageType, ImageType].New(maurer,
# Level=watershed_level_normalized, MarkWatershedLine=False)
# mask = itk.MaskImageFilter[ImageType, ImageType, ImageType].New(watershed, image)

stats = itk.BinaryImageToStatisticsLabelMapFilter[ImageType, ImageType, LabelMapType].New(image,image)
size = itk.ShapeOpeningLabelMapFilter[LabelMapType].New(stats,Attribute="PhysicalSize",
Lambda=physical_size_lambda)
relabel = itk.StatisticsRelabelLabelMapFilter[LabelMapType].New(size, Attribute='Mean')
label_to_binary = itk.LabelMapToBinaryImageFilter.New(Input=relabel)

# Fill Holes
filler = itk.BinaryFillholeImageFilter.New(Input=label_to_binary)
filler.FullyConnectedOn()
# TEST, comparer
# comparer = itk.AbsoluteValueDifferenceImageFilter.New(Input1=label_to_binary, Input2=filler)
# output_comparer = "/tmp/compare.nrrd"
# itk.ImageFileWriter.New(Input=comparer, FileName=output_comparer).Update()

itk.ImageFileWriter.New(Input=filler, FileName=output_filename).Update()

from subprocess import Popen
testDriver = "/home/phc/tmp/IsotropicWaveletsTestDriver"
pin = Popen([testDriver, 'runViewImage', input_filename, 'input'])
pout = Popen([testDriver, 'runViewImage', output_filename, 'relabel'])
# pout_inverted = Popen([testDriver, 'runViewImage', output_comparer, 'compare fill holes'])
