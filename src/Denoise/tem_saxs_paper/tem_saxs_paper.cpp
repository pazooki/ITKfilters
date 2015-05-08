
#include "gmock/gmock.h"
#include "Denoise.h"
#include <memory>
#include <QuickView.h>
using namespace testing;
using namespace std;

TEST(denoise, pectin1M1045){
    const string img{"./fixturesTemSaxsPaper/pectin_1_ice_Montage_1045_12K_8_bit.tif"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    // Visualize
    QuickView viewer;
    Denoise::InputImageType * pimg = r;
    viewer.AddImage(pimg, 1, img);

    auto outSmart = denoise->AnisotropicFilter(r, 5, 0.044, 2);
    Denoise::RealImageType* out = outSmart.GetPointer();
    viewer.AddImage( out, 1 , "Anisotropic Filter");

    auto outSmart2 = denoise->AnisotropicFilter(outSmart, 5, 0.044, 5, Denoise::AnisotropicFilterID::CURVATURE);
    Denoise::RealImageType* out2 = outSmart2.GetPointer();
    viewer.AddImage( out2, 1 , "Anisotropic Filterx2");

    Denoise::RealImageType* out3 = denoise->AnisotropicFilter(outSmart2, 10, 0.044, 2);
    viewer.AddImage( out3, 1 , "Anisotropic Filterx3");
    viewer.Visualize();
}
