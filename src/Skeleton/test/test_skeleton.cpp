#include "gmock/gmock.h"
#include <memory>
#include "../Denoise/test/prog_options_test.h"

bool VFLAG;
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    auto option_map = program_options(argc, argv);
    VFLAG = option_map["visualize"].as<bool>();
    return RUN_ALL_TESTS();
}

#include "Skeleton.h"
#include <QuickView.h>
#include "itkRandomImageSource.h"
using namespace testing;
using namespace std;

TEST(skeleton, works){
    const string img{"./fixtures/M1045_11_20.tiff"};
    cout << "Test works" << endl;
    // Visualize
    // QuickView viewer;
    // viewer.AddImage(r.GetPointer(), 1, img);
    // viewer.AddImage( outSmart->GetOutput(), 1 , "Anisotropic Filter");
    // auto outSmart2 = denoise->AnisotropicFilterCurvature(outSmart->GetOutput(), 5, 0.044, 5);
    // viewer.AddImage( outSmart2->GetOutput(), 1 , "Anisotropic Filterx2");
    // auto outSmart3 = denoise->AnisotropicFilterCurvature(outSmart2->GetOutput(), 10, 0.044, 2);
    // viewer.AddImage( outSmart3->GetOutput(), 1 , "Anisotropic Filterx3");
    // if (VFLAG) viewer.Visualize();
}

