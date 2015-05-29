
#include "gmock/gmock.h"
#include "Denoise.h"
#include <memory>
#include <QuickView.h>
#include "itkRandomImageSource.h"
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
TEST(composite, StructuringElement){
    typedef itk::FlatStructuringElement<2> FSL;
    typedef FSL::LType lt;
    // int objWidth = 2;
    // itk::SizeValueType
    // r1      = 5 + 2*objWidth;
    // auto r2      = 5 + 2*objWidth;
    FSL::RadiusType r;
    r[0]=5;
    r[1]=5;
    FSL el;
    el.SetRadius(r);
    el.SetDecomposable(1);
    float l = 3;
    float h = 2;
    // lt p10;
    // p10[0] = 0;
    // p10[1] = 0;
    // el.AddLine(p10);
    lt p11;
    p11[0] = l;
    p11[1] = 0;
    el.AddLine(p11);

    lt p20;
    p20[0] = 0;
    p20[1] = h;
    el.AddLine(p20);
    // lt p21;
    // p21[0] = l;
    // p21[1] = h;
    // el.AddLine(p21);
    el.ComputeBufferFromLines();

    typedef itk::Image< unsigned char, 2 >    ImageType;
    typedef itk::RandomImageSource<ImageType> RandomType;
    ImageType::SizeType size;
    size[0] = 20;
    size[1] = 30;
    auto image = RandomType::New();
    image->SetSize(size);
    image->Update();

    typedef itk::GrayscaleMorphologicalOpeningImageFilter<
        ImageType , ImageType, FSL >  OpeningFilterType;
    auto opening = OpeningFilterType::New();
    opening->SetKernel(el);
    opening->SetInput(image->GetOutput());
    opening->Update();
    QuickView viewer;
    viewer.AddImage(image->GetOutput());
    viewer.AddImage(opening->GetOutput());

    viewer.Visualize();
    // auto p2 = {{0,h},{l,h}};
    // el.Addline(p2);
    // auto f1 = {{0,0},{0,h}};
    // el.Addline(f1);
    // auto f2 = {{0,l},{l,h}};
    // el.Addline(f2);
    // auto lines = el.GetLines();
    // std::cout << lines[0][0] << " " << lines[0][1] << "; " << lines[1][0] <<" " <<lines[1][1] <<";" <<lines[2][0] << " " << lines[2][1] << ";" << lines[3][0] << " " << lines[3][1] << ";" << lines[4][0] << " " << lines[4][1] << endl;

    // for(int j=0; j<this->Size(); j++ )
    // {
    //     if( this->GetElement( j ) )
    //     {
    //         image->SetPixel( centerIdx+this->GetOffset( j ), foreground );
    //     }
    //     else
    //     {
    //         image->SetPixel( centerIdx+this->GetOffset( j ), background );
    //     }
    // }
    //
    // return image;

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
