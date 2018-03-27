import itk
import sys
import os
if len(sys.argv) < 3:
    print("Usage: " + sys.argv[0] + "inputImage outputFOLDER [maximum = 255]")
    sys.exit(1)
print("invert_intensity %s" % sys.argv[1])

input_filename = sys.argv[1]
basename = os.path.basename(input_filename)
filename, file_extension = os.path.splitext(basename)
output_folder = sys.argv[2]
maximum = 255;
if len(sys.argv) == 4:
    maximum = sys.argv[3]

output_filename = os.path.join(output_folder, filename +
                               "_INV" + file_extension)
# PixelType = itk.F
# Dimension = 3
# ImageType = itk.Image[PixelType, Dimension]
# reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader = itk.ImageFileReader.New(FileName=input_filename)
reader.Update()
# image = reader.GetOutput()
inverter = itk.InvertIntensityImageFilter.New(reader.GetOutput())
# If image is float, change int to float here. TODO: anyway to cast to the right type automatically?
inverter.SetMaximum(int(maximum))
inverter.Update()

itk.ImageFileWriter.New(Input=inverter, FileName=output_filename).Update()
