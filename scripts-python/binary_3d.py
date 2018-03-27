import itk
import sys
import os
if len(sys.argv) < 4:
  print("Usage: " + sys.argv[0] + "inputImage outputFOLDER otsu|0.0 invert")
  sys.exit(1)
print("binary_3d %s" % sys.argv[1])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename, file_extension = os.path.splitext(basename)
output_folder = sys.argv[2]
binary_method = sys.argv[3]
invert = "no"
if len(sys.argv) == 5:
  invert = "yes"
output_extension = "nrrd"
output_filename = os.path.join(output_folder, filename +
                               "_binary_" + ("_inverted_" if invert == "yes" else "") + binary_method + "." + output_extension)
PixelType = itk.F
Dimension = 3
ImageType = itk.Image[PixelType, Dimension]
reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()
binarizer = itk.OtsuThresholdImageFilter.New(image) if binary_method == "otsu" else itk.BinaryThresholdImageFilter.New(image)
if binary_method != "otsu":
  min_max_calculator = itk.MinimumMaximumImageCalculator.New(image)
  min_max_calculator.Compute()
  max_value = min_max_calculator.GetMaximum()
  min_value = min_max_calculator.GetMinimum()
  lower_threshold = (max_value - min_value)*float(binary_method)
  # binarizer.SetUpperThreshold(max_value);
  binarizer.SetLowerThreshold(lower_threshold);


if invert == "yes":
   inverter = itk.InvertIntensityImageFilter.New(binarizer)
   itk.ImageFileWriter.New(Input=inverter, FileName=output_filename).Update()
else:
  itk.ImageFileWriter.New(Input=binarizer, FileName=output_filename).Update()

from subprocess import Popen
testDriver = "/home/phc/tmp/IsotropicWaveletsTestDriver"
pin = Popen([testDriver, 'runViewImage', input_filename, 'input'])
pout = Popen([testDriver, 'runViewImage', output_filename, 'binary'])
# pout_inverted = Popen([testDriver, 'runViewImage', output_filename_inverted, 'otsu-inverted'])
