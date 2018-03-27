import prox_tv as ptv
import sys
import os
import itk
import numpy as np
import pywt

if len(sys.argv) != 4:
  print("Usage: " + sys.argv[0] + "inputImage outputFOLDER lambda")
  sys.exit(1)
print("TV denoising %s" % sys.argv[1])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename = os.path.splitext(basename)[0]
output_folder = sys.argv[2]
output_extension = "nrrd"
lam = float(sys.argv[3])
output_filename = os.path.join(output_folder, filename + "_tv" + str(lam) + "." + output_extension)
# PixelType = itk.F
# PixelType = itk.UC
# Dimension = 3
# ImageType = itk.Image[PixelType, Dimension]
# reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader = itk.ImageFileReader.New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()

image_array = itk.GetArrayViewFromImage(image)
norm = 1
dimensions_to_penalyze = [1,2,3]
tv_array = ptv.tvgen(image_array, np.array([lam,lam,lam]), dimensions_to_penalyze, np.array([norm,norm,norm]))

tv_array = np.ascontiguousarray(tv_array, dtype=image_array.dtype)
# modifiedImage = itk.PyBuffer[ImageType].GetImageViewFromArray(tv_array)
modifiedImage = itk.GetImageViewFromArray(tv_array)
itk.ImageFileWriter.New(Input=modifiedImage, FileName=output_filename).Update()

from subprocess import Popen
testDriver = "/home/phc/tmp/IsotropicWaveletsTestDriver"
pin = Popen([testDriver, 'runViewImage', input_filename, 'input'])
pout = Popen([testDriver, 'runViewImage', output_filename, 'TV'])
