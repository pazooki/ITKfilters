import itk
import sys
import os
if len(sys.argv) != 6:
  print("Usage: " + sys.argv[0] + "inputImage outputFOLDER majority neighborhood_radius iterations")
  sys.exit(1)
  print("Recommended: majority = 3, radius = 1, iterations = 10000")
print("binary_denoise_3d_fillholes_iterative %s" % sys.argv[1])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename = os.path.splitext(basename)[0]
output_folder = sys.argv[2]
# Majority is the number of pixels in the neighborhood of an OFF pixel, to turn it into ON.
# By default majority = 1, this means that an off pixel will be turned on if in the neighborhood (set by radius) there are at least 50% + 1 pixels ON.
# If radius = 1,1,1, neighborhood size will be 3x3 = 9 pixels.
# if 5 pixels around an OFF pixel are ON, then it will be switched.
majority = sys.argv[3]
radius = sys.argv[4]
iterations = sys.argv[5]
output_extension = "nrrd"
output_filename = os.path.join(output_folder, filename +
                               "_Majority" + majority +
                               "_Radius" + radius +
                               "_Iterations" + iterations +
                               "." + output_extension)
radius = float(radius)

PixelType = itk.UC
Dimension = 3
ImageType = itk.Image[PixelType, Dimension]
LabelObjectType = itk.StatisticsLabelObject[itk.UL, Dimension]
LabelMapType = itk.LabelMap[LabelObjectType]
reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()

# Fill Holes
filler = itk.VotingBinaryIterativeHoleFillingImageFilter.New(Input=image, MaximumNumberOfIterations=int(iterations), MajorityThreshold=int(majority))
radius_array = filler.GetRadius()
radius_array.Fill(int(radius))
filler.SetRadius(radius_array)
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
