
import itk
input_filename = "/home/phc/repository_local/ITKfilters/src/fixtures/collagen_32x32x16.tif"
output_filename = "/home/phc/repository_local/ITKfilters/scripts-python/results/collagen_32x32x16_otsu.tif"
PixelType = itk.F
Dimension = 3
ImageType = itk.Image[PixelType, Dimension]
reader = itk.ImageFileReader[ImageType].New(FileName=input_filename)
reader.Update()
image = reader.GetOutput()
binarizer = itk.OtsuThresholdImageFilter.New(image)

itk.ImageFileWriter.New(Input=binarizer, FileName=output_filename).Update()


#Python sucks
