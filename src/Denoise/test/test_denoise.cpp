
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
    denoise->AnisotropicFilterCurvature(r);
}

TEST(denoise, anisotropicFilterInPectinSubSet){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto outSmart = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 2);
    // // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage( outSmart->GetOutput(), 1 , "Anisotropic Filter");
    auto outSmart2 = denoise->AnisotropicFilterCurvature(outSmart->GetOutput(), 5, 0.044, 5);
    viewer.AddImage( outSmart2->GetOutput(), 1 , "Anisotropic Filterx2");
    auto outSmart3 = denoise->AnisotropicFilterCurvature(outSmart2->GetOutput(), 10, 0.044, 2);
    viewer.AddImage( outSmart3->GetOutput(), 1 , "Anisotropic Filterx3");
    viewer.Visualize();
}

TEST(denoise, RegionGrowthConnectedThreshold){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto d1 = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 5);
    auto rg = denoise->RegionGrowth(d1->GetOutput(), 0, 100);
    itk::Index<2> s1 = {{309,185}};
    rg->SetSeed(s1);
    rg->Update();
    auto rgo  = rg->GetOutput();
    // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(d1->GetOutput(), 1, "Denoised");
    viewer.AddImage(rgo, 1, "Rgrowth");
    viewer.Visualize();
}

TEST(denoise, write){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto d1 = denoise->AnisotropicFilterCurvature(r, 10, 0.044, 2);
    string outFile = "./testResults/M1045_11_20_denoised.tiff";
    denoise->Write(d1->GetOutput(), outFile);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(d1->GetOutput(), 1, "Denoised");
    viewer.Visualize();
}

TEST(morphological, open){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto mo10 = denoise->MorphologicalOpening(r, 10);
    auto mo5 = denoise->MorphologicalOpening(r, 5);
    auto mo20 = denoise->MorphologicalOpening(r, 20);
    auto mo3 = denoise->MorphologicalOpening(r, 3);
    // string outFile = "./testResults/M1045_11_20_morphOpen.tiff";
    // denoise->Write(mo, outFile);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(mo10->GetOutput(), 1, "Morphological Opening, r=10");
    viewer.AddImage(mo5->GetOutput(), 1, "Morphological Opening, r=5");
    viewer.AddImage(mo20->GetOutput(), 1, "Morphological Opening, r=20");
    viewer.AddImage(mo3->GetOutput(), 1, "Morphological Opening, r=3");
    auto d1 = denoise->AnisotropicFilterCurvature(mo3->GetOutput(), 5, 0.044, 5);
    viewer.AddImage(d1->GetOutput(), 1, "Morphological Opening + AnisotropicFilterCurvature, r=3");
    viewer.Visualize();
}
TEST(morphological, anisotropicPlusOpen){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto d1 = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 5);
    auto mo10 = denoise->MorphologicalOpening(d1->GetOutput(), 10);
    auto mo5 = denoise->MorphologicalOpening(d1->GetOutput(), 5);
    auto mo20 = denoise->MorphologicalOpening(d1->GetOutput(), 20);
    auto mo3 = denoise->MorphologicalOpening(d1->GetOutput(), 3);
    auto d2 = denoise->AnisotropicFilterCurvature(mo3->GetOutput(), 10, 0.044, 3);
    // string outFile = "./testResults/M1045_11_20_morphOpen.tiff";
    // denoise->Write(mo, outFile);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(d1->GetOutput(), 1, "Denoised");
    viewer.AddImage(mo10->GetOutput(), 1, "Morphological Opening, r=10");
    viewer.AddImage(mo5->GetOutput(), 1, "Morphological Opening, r=5");
    viewer.AddImage(mo20->GetOutput(), 1, "Morphological Opening, r=20");
    viewer.AddImage(mo3->GetOutput(), 1, "Morphological Opening, r=3");
    viewer.AddImage(d2->GetOutput(), 1, "Opening+Denoised, r=3");
    viewer.Visualize();
}
TEST(binary, Otsu){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryOtsu(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryOtsu");
    viewer.Visualize();
}
TEST(binary, Huang){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryHuang(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryHuang");
    viewer.Visualize();
}
TEST(binary, Yen){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryYen(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryYen");
    viewer.Visualize();
}
TEST(binary, Shanbhag){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryShanbhag(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryShanbhag");
    viewer.Visualize();
}
