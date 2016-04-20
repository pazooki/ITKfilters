#include "gtest/gtest.h"
#include <memory>
#include <string>
#include "prog_options_test.h"
#include "itkQuadratureFilterImageSource.h"
#include "itkImageFileReader.h"
#include "itkImage.h"

bool VFLAG;
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    auto option_map = program_options(argc, argv);
    VFLAG = option_map["visualize"].as<bool>();
    return RUN_ALL_TESTS();
}
using namespace std;
using namespace itk;

TEST(quadrature, InCollagen){
    const string img_file{"./fixtures/collagen_101x99x20.tiff"};
    const unsigned int dimension = 3;
    using PixelType = unsigned int;
    using ImageType = itk::Image<PixelType, dimension>;
    using ReaderType = itk::ImageFileReader<ImageType>;
    auto reader = ReaderType::New();
    reader->SetFileName(img_file);
    auto img = reader->GetOutput();
}
