import itk
import sys
import os
if len(sys.argv) < 4:
  print("Usage: " + sys.argv[0] + "outputImage geom.x geom.y geom.z inputImage1 inputImage2 ...")
  sys.exit(1)
print("montage_tiles %s" % sys.argv[1])

output_filename = sys.argv[1]
geom_x = sys.argv[2]
geom_y = sys.argv[3]
geom_z = sys.argv[4]
required_images=int(geom_x) * int(geom_y) * int(geom_z)
input_filename = sys.argv[5]
reader = itk.ImageFileReader.New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()

montager = itk.TileImageFilter.New(Input=image)
layout = montager.GetLayout()
layout[0] = int(geom_x)
layout[1] = int(geom_y)
layout[2] = int(geom_z)
montager.SetLayout(layout)

# More images:
print(required_images)
for im in range(1,required_images,1):
  input_filename = sys.argv[im + 5]
  reader = itk.ImageFileReader.New(FileName=input_filename)
  reader.Update()
  montager.SetInput(im, reader.GetOutput())


itk.ImageFileWriter.New(Input=montager, FileName=output_filename).Update()

from subprocess import Popen
testDriver = "/home/phc/tmp/IsotropicWaveletsTestDriver"
# pin = Popen([testDriver, 'runViewImage', input_filename, 'input'])
pout = Popen([testDriver, 'runViewImage', output_filename, 'montage'])
# pout_inverted = Popen([testDriver, 'runViewImage', output_comparer, 'compare fill holes'])
