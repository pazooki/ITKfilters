import sys
import os
import itk
import numpy as np

if len(sys.argv) < 5 or len(sys.argv) > 6:
  print("Usage: " + sys.argv[0] + "inputImage outputFOLDER numberOfIterations conductance timeStep  ")
  sys.exit(1)
print("Anisotropic Denoising %s" % sys.argv[1])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename = os.path.splitext(basename)[0]
output_folder = sys.argv[2]
output_extension = "nrrd"
iters = int(sys.argv[3])
conductance = float(sys.argv[4])
time_step = 1.0/(2.0**4) # 0.0625
if len(sys.argv) == 5:
  time_step = float(sys.argv[4])
output_filename = os.path.join(output_folder, filename +
                               "_AnisDenoise_c" + str(conductance) +
                               "_N" + str(iters) +
                               "_t" + str(time_step) +
                               "." + output_extension)
PixelType = itk.F
Dimension = 3
ImageType = itk.Image[PixelType, Dimension]
reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader.Update()

denoiser = itk.CurvatureAnisotropicDiffusionImageFilter.New(Input=reader.GetOutput())
denoiser.Update()

itk.ImageFileWriter.New(Input=denoiser.GetOutput(), FileName=output_filename).Update()

from subprocess import Popen
testDriver = "/home/phc/tmp/IsotropicWaveletsTestDriver"
pin = Popen([testDriver, 'runViewImage', input_filename, 'input'])
pout = Popen([testDriver, 'runViewImage', output_filename, 'denoised'])
