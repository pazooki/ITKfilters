#include "Denoise.h"
#include "gmock/gmock.h"
#include <memory>
#include <QuickView.h>
#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageDuplicator.h>
#include <algorithm>
#include <itkPasteImageFilter.h>
#include <itkExtractImageFilter.h>
using namespace testing;
using namespace std;
TEST(binary, pectin1M1045BinaryGlobal) {
    const string img{"./resultsTemSaxsPaper/pectin1_1045_homogeneous.tif"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);

    auto d1 = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 5);
    auto mo3 = denoise->MorphologicalOpening(d1->GetOutput(), 3);
    auto d2 = denoise->AnisotropicFilterCurvature(mo3->GetOutput(), 10, 0.044, 3);
    auto bOtsu = denoise->BinaryOtsu(d2->GetOutput());
    auto moBinary = denoise->MorphologicalOpening(bOtsu->GetOutput(), 4);
    auto mcBinary = denoise->MorphologicalClosing(moBinary->GetOutput(), 4);
    // auto bHuang = denoise->BinaryHuang(d2);
    // auto bYen = denoise->BinaryHuang(d2);
    // auto bShanbhag = denoise->BinaryShanbhag(d2);
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);
    viewer.AddImage(d1->GetOutput(), 1, "Denoised");
    viewer.AddImage(mo3->GetOutput(), 1, "Morphological Opening, r=3");
    viewer.AddImage(d2->GetOutput(), 1, "Opening+Denoised, r=3");
    viewer.AddImage(bOtsu->GetOutput(), 1, "BinaryOtsu Denoised+Opening+Denoised");
    viewer.AddImage(moBinary->GetOutput(), 1, "Binary Opened r=4");
    // viewer.AddImage(bHuang.GetPointer(), 1, "BinaryHuang Denoised+Opening+Denoised");
    // viewer.AddImage(bYen.GetPointer(), 1, "BinaryYen Denoised+Opening+Denoised");
    // viewer.AddImage(bShanbhag.GetPointer(), 1, "BinaryShanbhag Denoised+Opening+Denoised");

    string outFile = "./resultsTemSaxsPaper/pectin1_1045_homogeneous_binarized.tif";
    denoise->Write(bOtsu->GetOutput(),outFile);
    outFile = "./resultsTemSaxsPaper/pectin1_1045_homogeneous_binarized_opened_closed.tif";
    denoise->Write(mcBinary->GetOutput(),outFile);
    viewer.Visualize();
}

TEST(binary, pectin1M1045BinaryRegions) {
    const string img{"./resultsTemSaxsPaper/pectin1_1045_homogeneous.tif"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    // Set Region
    Denoise::InputImageType::RegionType region;
    Denoise::InputImageType::SizeType blockSize;
    auto imgSize = denoise->inputImg_->GetLargestPossibleRegion().GetSize();
    int Nimages    = 5;
    blockSize[0]   = imgSize[0]/Nimages;
    blockSize[1]   = imgSize[1]/Nimages;
    region.SetSize(blockSize);
    Denoise::InputImageType::IndexType indexStart;
    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto cloneFilter = CastFilterType::New();
    cloneFilter->SetInput(r);
    cloneFilter->Update();
    Denoise::RealImageType::Pointer clone = cloneFilter->GetOutput();
    typedef itk::PasteImageFilter <Denoise::RealImageType, Denoise::RealImageType >  PasteImageFilterType;
    typedef itk::ExtractImageFilter< Denoise::InputImageType, Denoise::RealImageType > ExtractFilterType;
    for (int i = 0; i < Nimages ; i++){
        for (int j = 0; j < Nimages ; j++){

            indexStart[0] = blockSize[0] * i;
            indexStart[1] = blockSize[1] * j;
            region.SetIndex(indexStart);
            auto crop = ExtractFilterType::New();
            crop->SetExtractionRegion(region);
            crop->SetInput(r);
            crop->Update();
            auto d1 = denoise->AnisotropicFilterCurvature(crop->GetOutput(), 5, 0.044, 5);
            auto mo3 = denoise->MorphologicalOpening(d1->GetOutput(), 3);
            auto d2 = denoise->AnisotropicFilterCurvature(mo3->GetOutput(), 10, 0.044, 3);
            auto bHuang = denoise->BinaryOtsu(d2->GetOutput());
            auto moB = denoise->MorphologicalOpening(bHuang->GetOutput(), 4);
            auto mocB = denoise->MorphologicalClosing(moB->GetOutput(), 2);
            mocB->Update();
            // cout << "beforePaste " << i << " " << j << endl;
            // auto s = mocB->GetOutput()->GetLargestPossibleRegion().GetSize();
            // cout << s[0] << s[1] << endl;
            // s = clone->GetLargestPossibleRegion().GetSize();
            // cout << s[0] << s[1] << endl;
            PasteImageFilterType::Pointer pasteFilter = PasteImageFilterType::New();
            pasteFilter->SetSourceImage(mocB->GetOutput());
            pasteFilter->SetDestinationImage(clone);
            pasteFilter->SetSourceRegion(mocB->GetOutput()->GetLargestPossibleRegion());
            pasteFilter->SetDestinationIndex(indexStart);
            pasteFilter->Update();
            clone = pasteFilter->GetOutput();
            clone->DisconnectPipeline();
            // cout << "afterPaste " << i << " " << j <<endl;
            // s = mocB->GetOutput()->GetLargestPossibleRegion().GetSize();
            // cout << s[0] << s[1] << endl;
            // s = clone->GetLargestPossibleRegion().GetSize();
            // cout << s[0] << s[1] << endl;
        }
    }
    string outFile = "./resultsTemSaxsPaper/pectin1_1045_binary_local.tif";
    denoise->Write(clone,outFile);
}

TEST(homogenous, pectin1M1045){
    const string img{"./resultsTemSaxsPaper/pectin_1_ice_Montage_1045_12K_8_bit.tif"};
    auto denoise = make_shared<Denoise>();
    auto r       = denoise->Read(img);
    Denoise::InputImageType::RegionType region;
    Denoise::InputImageType::SizeType blockSize;
    auto imgSize = denoise->inputImg_->GetLargestPossibleRegion().GetSize();
    int Nimages    = 5;
    blockSize[0]   = imgSize[0]/Nimages;
    blockSize[1]   = imgSize[1]/Nimages;
    region.SetSize(blockSize);

    Denoise::InputImageType::IndexType indexStart;
    // The montage is made by 5x5 images, iterate over all of them.
    vector<double> intensityRegion;
    for (int i = 0; i < Nimages ; i++){
        for (int j = 0; j < Nimages ; j++){

            indexStart[0] = blockSize[0] * i;
            indexStart[1] = blockSize[1] * j;

            region.SetIndex(indexStart);

            itk::ImageRegionConstIterator<Denoise::InputImageType> cIter(r, region);
            // Get the Mean Intensity of the Region (normalized by the size)
            int regionIntensity = 0;
            unsigned char val;
            while(!cIter.IsAtEnd())
            {
                val = cIter.Get();
                regionIntensity += val;
                ++cIter;
            }
            double regionIntensityD = regionIntensity / (blockSize[0] * blockSize[1]);
            intensityRegion.push_back(regionIntensityD);
        }
    }

    // Mean value of Intensity:
    double mean = 0.0;
    for (auto v : intensityRegion){
       mean += v;
    }
    mean /= intensityRegion.size();
    cout << mean << endl;

    vector<double> correctFactor = intensityRegion;
    for (auto &v : correctFactor){
        v = v/mean;
        // v = v / *min_element(intensityRegion.begin(), intensityRegion.end());
    }

    // Create a duplicate of input image and modify the intensity of the region blocks:
    typedef itk::ImageDuplicator< Denoise::InputImageType > DuplicatorType;
    auto duplicator = DuplicatorType::New();
    duplicator->SetInputImage(denoise->inputImg_);
    duplicator->Update();
    Denoise::InputImageType::Pointer clone = duplicator->GetModifiableOutput();

    for (int i = 0; i < Nimages ; i++){
        for (int j = 0; j < Nimages ; j++){

            indexStart[0] = blockSize[0] * i;
            indexStart[1] = blockSize[1] * j;

            region.SetIndex(indexStart);

            itk::ImageRegionIterator<Denoise::InputImageType> iter(clone, region);
            // Get the Mean Intensity of the Region (normalized by the size)
            unsigned short val;
            unsigned short nval;
            while(!iter.IsAtEnd())
            {
                val  = iter.Get();
                nval = val / correctFactor[i*Nimages + j];
                // cout << nval << endl;
                iter.Set(nval);
                ++iter;
            }
        }
    }
    string outFile = "./resultsTemSaxsPaper/pectin1_1045_homogeneous.tif";
    denoise->Write(clone,outFile);

}

TEST(homogeneous, pectin1M1045PlusDenoise) {
    const string img{"./resultsTemSaxsPaper/pectin1_1045_homogeneous.tif"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);

    auto out = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 2);
    viewer.AddImage( out->GetOutput(), 1 , "Anisotropic Filter");

    auto out2 = denoise->AnisotropicFilterCurvature(out->GetOutput(), 5, 0.044, 5);
    viewer.AddImage( out2->GetOutput(), 1 , "Anisotropic Filterx2");

    auto out3 = denoise->AnisotropicFilterCurvature(out2->GetOutput(), 10, 0.044, 2);
    viewer.AddImage( out3->GetOutput(), 1 , "Anisotropic Filterx3");
    string outFile = "./resultsTemSaxsPaper/pectin1_1045_homogeneous_denoised.tif";
    denoise->Write(out3->GetOutput(),outFile);
    // viewer.Visualize();
}

TEST(denoise, pectin1M1045){
    const string img{"./resultsTemSaxsPaper/pectin_1_ice_Montage_1045_12K_8_bit.tif"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);

    auto out = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 2);
    viewer.AddImage( out->GetOutput(), 1 , "Anisotropic Filter");

    auto out2 = denoise->AnisotropicFilterCurvature(out->GetOutput(), 5, 0.044, 5);
    viewer.AddImage( out2->GetOutput(), 1 , "Anisotropic Filterx2");

    auto out3 = denoise->AnisotropicFilterCurvature(out2->GetOutput(), 10, 0.044, 2);
    viewer.AddImage( out3->GetOutput(), 1 , "Anisotropic Filterx3");
    string outFile = "./resultsTemSaxsPaper/pectin1_1045_denoised.tif";
    denoise->Write(out3->GetOutput(),outFile);
    viewer.Visualize();
}

TEST(denoise, carrageenanK832){
    const string img{"./resultsTemSaxsPaper/Montage_832.tif"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);

    auto out = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 2);
    viewer.AddImage( out->GetOutput(), 1 , "Anisotropic Filter");

    auto out2 = denoise->AnisotropicFilterCurvature(out->GetOutput(), 5, 0.044, 5);
    viewer.AddImage( out2->GetOutput(), 1 , "Anisotropic Filterx2");

    auto out3 = denoise->AnisotropicFilterCurvature(out2->GetOutput(), 10, 0.044, 2);
    viewer.AddImage( out3->GetOutput(), 1 , "Anisotropic Filterx3");
    // viewer.Visualize();
    string outFile = "./resultsTemSaxsPaper/CaKMontage832_denoised.tif";
    denoise->Write(out3->GetOutput(),outFile);
}

TEST(denoise, carrageenanNa851){
    const string img{"./resultsTemSaxsPaper/Montage_851.tif"};
    auto denoise = make_shared<Denoise>() ;
    auto r = denoise->Read(img);
    // Visualize
    QuickView viewer;
    viewer.AddImage(r.GetPointer(), 1, img);

    auto out = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 2);
    viewer.AddImage( out->GetOutput(), 1 , "Anisotropic Filter");

    auto out2 = denoise->AnisotropicFilterCurvature(out->GetOutput(), 5, 0.044, 5);
    viewer.AddImage( out2->GetOutput(), 1 , "Anisotropic Filterx2");

    auto out3 = denoise->AnisotropicFilterCurvature(out2->GetOutput(), 10, 0.044, 2);
    viewer.AddImage( out3->GetOutput(), 1 , "Anisotropic Filterx3");
    // viewer.Visualize();
    string outFile = "./resultsTemSaxsPaper/CaNaMontage851_denoised.tif";
    denoise->Write(out3->GetOutput(),outFile);
}

// TEST(binary, pectin1M1045BinaryOpenCloseOpen){
//     const string img{"./resultsTemSaxsPaper/pectin1_1045_homogeneous.tif"};
//     auto denoise = make_shared<Denoise>() ;
//     auto r = denoise->Read(img);
//     auto d1 = denoise->AnisotropicFilterCurvature(r, 5, 0.044, 5);
//     auto mo3 = denoise->MorphologicalOpening(d1, 3);
//     auto d2 = denoise->AnisotropicFilterCurvature(mo3, 10, 0.044, 3);
//     auto bHuang = denoise->BinaryHuang(d2);
//     auto moB = denoise->MorphologicalOpening(bHuang, 4);
//     auto mocB = denoise->MorphologicalClosing(moB, 2);
//     auto mcB = denoise->MorphologicalClosing(bHuang, 1);
//     auto mcoB = denoise->MorphologicalOpening(mcB, 4);
//     // Visualize
//     QuickView viewer;
//     viewer.AddImage(r.GetPointer(), 1, img);
//     viewer.AddImage(d1.GetPointer(), 1, "Denoised");
//     viewer.AddImage(mo3.GetPointer(), 1, "Morphological Opening, r=3");
//     viewer.AddImage(d2.GetPointer(), 1, "Opening+Denoised, r=3");
//     viewer.AddImage(bHuang.GetPointer(), 1, "BinaryHuang Denoised+Opening+Denoised");
//     viewer.AddImage(moB.GetPointer(), 1, "Huang ... + Opening");
//     viewer.AddImage(mcB.GetPointer(), 1, "Huang ... + Closing");
//     viewer.AddImage(mocB.GetPointer(), 1, "Huang ... + Opening + Closing");
//     viewer.AddImage(mcoB.GetPointer(), 1, "Huang ... + Closing + Opening");
//     viewer.Visualize();
//     string outFile = "./resultsTemSaxsPaper/pectin1_1045_binaryHuangCloseOpen.tif";
//     denoise->Write(mcoB,outFile);
//     // viewer.Visualize();
// }