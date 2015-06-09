#include "gmock/gmock.h"
#include "itkFlatStructuringElementWithImageBridge.h"
#include <memory>
#include <QuickView.h>
#include "itkRandomImageSource.h"
#include <itkImageFileReader.h>
using namespace testing;
using namespace std;

TEST(ImageBridge, constructor){

    const static unsigned int  Dimension = 2;
    typedef unsigned int       PixelType;
    typedef itk::Image< PixelType, Dimension> ImageType;
    typedef itk::RandomImageSource<ImageType> RandomType;
    ImageType::SizeType size;
    size[0] = 5;
    size[1] = 5;
    auto rimage = RandomType::New();
    rimage->SetSize(size);
    rimage->Update();

    typedef itk::ImageFileReader< ImageType > ReaderType;
    auto reader = ReaderType::New();

    string inputName = "./fixtures/cyld3.png";
    reader->SetFileName( inputName );
    try {
       reader->Update();
    } catch( itk::ExceptionObject & excp ) {
        std::cerr << "Problem encountered while reading image file : " <<
            inputName << std::endl;
        throw;
    }

    typedef itk::FlatStructuringElementWithImageBridge<Dimension,ImageType> FlatBridgeType;
    auto f = FlatBridgeType();
    auto f2 = FlatBridgeType::FromImage(reader->GetOutput());

    QuickView viewer;
    viewer.AddImage(reader->GetOutput());
    viewer.AddImage(f2.GetImage().GetPointer());

    viewer.Visualize();
}
