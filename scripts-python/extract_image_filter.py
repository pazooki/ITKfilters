import itk
import sys
import os
if len(sys.argv) != 9:
  print("Usage: " + sys.argv[0] + "inputImage outputImage index.x index.y index.z size.x size.y size.z")
  sys.exit(1)
print("extract_image_filter %s" % sys.argv[1])
print("output: %s" % sys.argv[2])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename = os.path.splitext(basename)[0]
output_filename = sys.argv[2]
index_x = sys.argv[3]
index_y = sys.argv[4]
index_z = sys.argv[5]
size_x = sys.argv[6]
size_y = sys.argv[7]
size_z = sys.argv[8]

reader = itk.ImageFileReader.New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()

cropper = itk.ExtractImageFilter.New(Input=image)
cropper.SetDirectionCollapseToIdentity()
extraction_region = cropper.GetExtractionRegion()
size = extraction_region.GetSize()
size[0] = int(size_x)
size[1] = int(size_y)
size[2] = int(size_z)
index = extraction_region.GetIndex()
index[0] = int(index_x)
index[1] = int(index_y)
index[2] = int(index_z)
extraction_region.SetSize(size)
extraction_region.SetIndex(index)
cropper.SetExtractionRegion(extraction_region)

itk.ImageFileWriter.New(Input=cropper, FileName=output_filename).Update()

# from subprocess import Popen
# testDriver = "/home/phc/tmp/IsotropicWaveletsTestDriver"
# pin = Popen([testDriver, 'runViewImage', input_filename, 'input'])
# pout = Popen([testDriver, 'runViewImage', output_filename, 'cropped'])
# pout_inverted = Popen([testDriver, 'runViewImage', output_comparer, 'compare fill holes'])
