/* #include "gmock/gmock.h" */
#include "gtest/gtest.h"
#include <memory>
#include "prog_options_test.h"

bool VFLAG;
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    auto option_map = program_options(argc, argv);
    VFLAG = option_map["visualize"].as<bool>();
    return RUN_ALL_TESTS();
}

#include "Denoise.h"
#include <QuickView.h>
#include "itkRandomImageSource.h"
#include <itkImageFileWriter.h>
using namespace testing;
using namespace std;
TEST(fileIO, readTiff){
    const string img{"./fixtures/imgTiny.tiff"};
    auto denoise = make_shared<Denoise>() ;
    denoise->Read(img);
}
TEST(fileIO, readPng){
    const string img{"./fixtures/noiseRandom200x200.png"};
    auto denoise = make_shared<Denoise>() ;
    denoise->Read(img);
}

TEST(fileIO, write){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto d1 = denoise->AnisotropicFilterCurvature(r, 10, 0.044, 2);
    string outFile = "./testResults/M1045_11_20_denoised.tiff";
    denoise->Write(d1->GetOutput(), outFile);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(d1->GetOutput(), 1, "Denoised");
    if (VFLAG) viewer.Visualize();

}
TEST(anisotropic, InNoise){
    const string img{"./fixtures/noiseRandom200x200.png"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    denoise->AnisotropicFilterCurvature(r);
}

TEST(anisotropic, InPectinSubSet){
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
    if (VFLAG) viewer.Visualize();
}

TEST(regionGrowth, connectedThreshold){
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
    if (VFLAG) viewer.Visualize();
}


TEST(morphological, opening){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto mo10 = denoise->MorphologicalOpening(r, 10);
    auto mo5 = denoise->MorphologicalOpening(r, 5);
    auto mo20 = denoise->MorphologicalOpening(r, 20);
    auto mo3 = denoise->MorphologicalOpening(r, 3);

    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(mo10->GetOutput(), 1, "Morphological Opening, r=10");
    viewer.AddImage(mo5->GetOutput(), 1, "Morphological Opening, r=5");
    viewer.AddImage(mo20->GetOutput(), 1, "Morphological Opening, r=20");
    viewer.AddImage(mo3->GetOutput(), 1, "Morphological Opening, r=3");
    auto d1 = denoise->AnisotropicFilterCurvature(mo3->GetOutput(), 5, 0.044, 5);
    if (VFLAG) viewer.Visualize();
}
TEST(morphological, closing){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto mc10 = denoise->MorphologicalClosing(r, 10);
    auto mc5 = denoise->MorphologicalClosing(r, 5);
    auto mc20 = denoise->MorphologicalClosing(r, 20);
    auto mc3 = denoise->MorphologicalClosing(r, 3);

    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(mc10->GetOutput(), 1, "Morphological Opening, r=10");
    viewer.AddImage(mc5->GetOutput(), 1, "Morphological Opening, r=5");
    viewer.AddImage(mc20->GetOutput(), 1, "Morphological Opening, r=20");
    viewer.AddImage(mc3->GetOutput(), 1, "Morphological Opening, r=3");
    auto d1 = denoise->AnisotropicFilterCurvature(mc3->GetOutput(), 5, 0.044, 5);
    if (VFLAG) viewer.Visualize();
}

TEST(binary, threshold){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    unsigned int threshold = 90;
    auto b = denoise->BinaryThreshold(r, threshold);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryThreshold");
    if (VFLAG) viewer.Visualize();
    string outFile = "./testResults/pectin_threshold.tiff";
    denoise->Write(b->GetOutput(), outFile);
}
TEST(binary, Otsu){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryOtsu(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryOtsu");
    if (VFLAG) viewer.Visualize();
}
TEST(binary, Huang){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryHuang(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryHuang");
    if (VFLAG) viewer.Visualize();
}
TEST(binary, Yen){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryYen(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryYen");
    if (VFLAG) viewer.Visualize();
}
TEST(binary, Shanbhag){
    const string img{"./fixtures/M1045_11_20.tiff"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    auto b = denoise->BinaryShanbhag(r);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(b->GetOutput(), 1, "BinaryShanbhag");
    if (VFLAG) viewer.Visualize();
}
