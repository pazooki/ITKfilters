import itk
import sys
import os

if len(sys.argv) != 3:
  print("Usage: " + sys.argv[0] + "inputImage outputFOLDER")
  sys.exit(1)
print("Otsu binarizing %s" % sys.argv[1])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename = os.path.splitext(basename)[0]
output_folder = sys.argv[2]
output_extension = "nrrd"
output_filename = os.path.join(output_folder, filename + "_otsu." + output_extension)
PixelType = itk.F
Dimension = 3
ImageType = itk.Image[PixelType, Dimension]
reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()
binarizer = itk.OtsuThresholdImageFilter.New(image)
itk.ImageFileWriter.New(Input=binarizer, FileName=output_filename).Update()

output_filename_inverted = os.path.join(output_folder, filename + "_otsu_inverted." +  output_extension)
inverted = itk.InvertIntensityImageFilter.New(binarizer)
itk.ImageFileWriter.New(Input=inverted, FileName=output_filename_inverted).Update()

from subprocess import Popen
testDriver = "/home/phc/tmp/IsotropicWaveletsTestDriver"
pin = Popen([testDriver, 'runViewImage', input_filename, 'input'])
pout = Popen([testDriver, 'runViewImage', output_filename, 'otsu'])
pout_inverted = Popen([testDriver, 'runViewImage', output_filename_inverted, 'otsu-inverted'])
