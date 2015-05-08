
#include "gmock/gmock.h"
#include "Denoise.h"
#include <memory>
#include <QuickView.h>
using namespace testing;
using namespace std;

TEST(denoise, readTiff){
    const string img{"./fixtures/imgTiny.tiff"};
    auto denoise = make_shared<Denoise>() ;
    denoise->Read(img);
}
TEST(denoise, readPng){
    const string img{"./fixtures/noiseRandom200x200.png"};
    auto denoise = make_shared<Denoise>() ;
    denoise->Read(img);
}

TEST(denoise, anisotropicFilterInNoise){
    const string img{"./fixtures/noiseRandom200x200.png"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    denoise->AnisotropicFilter(r);
}

TEST(denoise, anisotropicFilterInPectinSubSet){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto outSmart = denoise->AnisotropicFilter(r, 5, 0.044, 2);
    // // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage( outSmart.GetPointer(), 1 , "Anisotropic Filter");
    auto outSmart2 = denoise->AnisotropicFilter(outSmart, 5, 0.044, 5, Denoise::AnisotropicFilterID::CURVATURE);
    viewer.AddImage( outSmart2.GetPointer(), 1 , "Anisotropic Filterx2");
    auto outSmart3 = denoise->AnisotropicFilter(outSmart2, 10, 0.044, 2);
    viewer.AddImage( outSmart3.GetPointer(), 1 , "Anisotropic Filterx3");
    viewer.Visualize();
}

TEST(denoise, RegionGrowthConnectedThreshold){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto d1 = denoise->AnisotropicFilter(r, 5, 0.044, 5, Denoise::AnisotropicFilterID::CURVATURE);
    auto rg = denoise->RegionGrowth(d1, 0, 100);
    itk::Index<2> s1 = {{309,185}};
    rg->SetSeed(s1);
    rg->Update();
    auto rgo  = rg->GetOutput();
    // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(d1.GetPointer(), 1, "Denoised");
    viewer.AddImage(rgo, 1, "Rgrowth");
    viewer.Visualize();
}
