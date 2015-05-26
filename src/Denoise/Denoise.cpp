#include "Denoise.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkCastImageFilter.h>
using namespace std;

Denoise::InputTypeP Denoise::Read(const string &inputName){
    // Reader
    typedef itk::ImageFileReader< InputImageType > ReaderType;
    auto reader = ReaderType::New();

    // SCIFIO type
    // typedef itk::SCIFIOImageIO  ImageIOType;
    // TIFF type
    // typedef itk::TIFFImageIO ImageIOType;
    // auto scifioIO = ImageIOType::New();
    // reader->SetImageIO(scifioIO);
    reader->SetFileName( inputName );

    try {
       reader->Update();
    } catch( itk::ExceptionObject & excp ) {
        std::cerr << "Problem encountered while reading image file : " <<
            inputName << std::endl;
        throw;
    }

    // typedef itk::CastImageFilter< InputImageType, InputImageType > castFilterType;
    // castFilterType::Pointer castFilter = castFilterType::New();
    // castFilter->SetInput(reader->GetOutput());
   // castFilter->Update();

    inputImg_ = reader->GetOutput();

    // QuickView viewer;
    // viewer.AddImage(inputImg_);
    // viewer.Visualize();
    // inputImg_->Register();
    return inputImg_;
}

void Denoise::Write(Denoise::RealImageType* input, std::string &imgName){

    typedef unsigned char                           OutputPixelType;
    typedef itk::Image< OutputPixelType, 2 > OutputImageType;
    typedef itk::RescaleIntensityImageFilter<Denoise::RealImageType, OutputImageType > RescaleFilter;
    auto rescaleFilter = RescaleFilter::New();
    rescaleFilter->SetInput(input);
    rescaleFilter->SetOutputMinimum( itk::NumericTraits< OutputPixelType >::min() );
    rescaleFilter->SetOutputMaximum( itk::NumericTraits< OutputPixelType >::max() );
    typedef  itk::ImageFileWriter< OutputImageType  > WriterType;
    auto writer = WriterType::New();
    writer->SetFileName(imgName);
    writer->SetInput(rescaleFilter->GetOutput());
    writer->Update();
    cout << "Output generated at " << imgName << endl;
}

void Denoise::Write(Denoise::InputImageType* input, std::string &imgName){

    typedef unsigned char                           OutputPixelType;
    typedef itk::Image< OutputPixelType, 2 > OutputImageType;
    typedef itk::RescaleIntensityImageFilter<Denoise::InputImageType, OutputImageType > RescaleFilter;
    auto rescaleFilter = RescaleFilter::New();
    rescaleFilter->SetInput(input);
    rescaleFilter->SetOutputMinimum( itk::NumericTraits< OutputPixelType >::min() );
    rescaleFilter->SetOutputMaximum( itk::NumericTraits< OutputPixelType >::max() );
    typedef  itk::ImageFileWriter< OutputImageType  > WriterType;
    auto writer = WriterType::New();
    writer->SetFileName(imgName);
    writer->SetInput(rescaleFilter->GetOutput());
    writer->Update();
    cout << "Output generated at " << imgName << endl;

}
void Denoise::Write(Denoise::RealTypeP &input, std::string &imgName){

    Write(input.GetPointer(), imgName);
}

void Denoise::Write(Denoise::InputTypeP &input, std::string &imgName){

    Write(input.GetPointer(), imgName);
}

Denoise::AnisotropicFilterGradientTypeP Denoise::AnisotropicFilterGradient(Denoise::RealTypeP img,
        const unsigned int numberOfIterations, const double timeStep, const double conductance)
{
    auto anisotropicFilter = AnisotropicFilterGradientType::New();
    anisotropicFilter->SetInput( img );
    // numberOfIterations = 5;
    // const double timeStep = 0.125;
    //  const double conductance = 5;
    anisotropicFilter->SetNumberOfIterations( numberOfIterations );
    anisotropicFilter->SetTimeStep( timeStep );
    anisotropicFilter->SetConductanceParameter( conductance );
//    anisotropicFilter->Update();
    return anisotropicFilter;
}

Denoise::AnisotropicFilterGradientTypeP Denoise::AnisotropicFilterGradient(Denoise::InputTypeP img,
        const unsigned int numberOfIterations, const double timeStep,
        const double conductance)
{

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    // InputImageType i = *img->Clone();
    filter->SetInput(img);
    filter->Update();
    return AnisotropicFilterGradient(filter->GetOutput(),
                numberOfIterations,timeStep, conductance);
}
Denoise::AnisotropicFilterCurvatureTypeP Denoise::AnisotropicFilterCurvature(Denoise::RealTypeP img,
        const unsigned int numberOfIterations, const double timeStep, const double conductance)
{
    auto anisotropicFilter = AnisotropicFilterCurvatureType::New();
    anisotropicFilter->SetInput( img );
    // numberOfIterations = 5;
    // const double timeStep = 0.125;
    //  const double conductance = 5;
    anisotropicFilter->SetNumberOfIterations( numberOfIterations );
    anisotropicFilter->SetTimeStep( timeStep );
    anisotropicFilter->SetConductanceParameter( conductance );
//    anisotropicFilter->Update();
    return anisotropicFilter;
}

Denoise::AnisotropicFilterCurvatureTypeP Denoise::AnisotropicFilterCurvature(Denoise::InputTypeP img,
        const unsigned int numberOfIterations, const double timeStep,
        const double conductance)
{

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    // InputImageType i = *img->Clone();
    filter->SetInput(img);
    filter->Update();
    return AnisotropicFilterCurvature(filter->GetOutput(),
                numberOfIterations,timeStep, conductance);
}

Denoise::ConnectedFilterTypeP Denoise::RegionGrowth(Denoise::RealTypeP img,
                    unsigned int lth, unsigned int hth)
{
    filterGrow_ = ConnectedFilterType::New();
    filterGrow_->SetInput(img);
    filterGrow_->SetLower( lth );
    filterGrow_->SetUpper( hth );
    filterGrow_->SetReplaceValue( 255 );
    return filterGrow_;
}

Denoise::ConnectedFilterTypeP Denoise::RegionGrowth(Denoise::InputTypeP img,
                    unsigned int lth, unsigned int hth)
{
    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    // InputImageType i = *img->Clone();
    filter->SetInput(img);
    filter->Update();
    return RegionGrowth(filter->GetOutput(),
                lth,hth);
}
Denoise::OpeningFilterTypeP Denoise::MorphologicalOpening(Denoise::RealTypeP img, int radius){

    // structuring element
    typedef itk::BinaryBallStructuringElement<
        unsigned char, 2  > StructuringElementType;
    // define the opening type
    typedef itk::GrayscaleMorphologicalOpeningImageFilter<
        RealImageType , RealImageType, StructuringElementType >  OpeningFilterType;
    // Create structuring element
    StructuringElementType  structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();
    // Opening filter
    auto opening = OpeningFilterType::New();
    opening->SetKernel(structuringElement);
    opening->SetInput(img);
//    opening->Update();
    return opening;
}
Denoise::OpeningFilterTypeP Denoise::MorphologicalOpening(Denoise::InputTypeP img, int radius){

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    filter->SetInput(img);
    filter->Update();
    return MorphologicalOpening(filter->GetOutput(), radius);
    // typedef itk::CastImageFilter< Denoise::RealImageType, Denoise::InputImageType > CastFilterTypeInverse;
    // auto filterInv = CastFilterTypeInverse::New();
    // // InputImageType i = *img->Clone();
    // filterInv->SetInput(mo);
    // filterInv->Update();
    // return filterInvGetOutput();

}

Denoise::ClosingFilterTypeP Denoise::MorphologicalClosing(Denoise::RealTypeP img, int radius){

    // structuring element
    typedef itk::BinaryBallStructuringElement<
        unsigned char, 2  > StructuringElementType;
    // define the opening type
    typedef itk::GrayscaleMorphologicalClosingImageFilter<
        RealImageType , RealImageType, StructuringElementType >  ClosingFilterType;
    // Create structuring element
    StructuringElementType  structuringElement;
    structuringElement.SetRadius(radius);
    structuringElement.CreateStructuringElement();
    // Closing filter
    auto opening = ClosingFilterType::New();
    opening->SetKernel(structuringElement);
    opening->SetInput(img);
//    opening->Update();
    return opening;
}
Denoise::ClosingFilterTypeP Denoise::MorphologicalClosing(Denoise::InputTypeP img, int radius){

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    filter->SetInput(img);
    filter->Update();
    return MorphologicalClosing(filter->GetOutput(), radius);
    // typedef itk::CastImageFilter< Denoise::RealImageType, Denoise::InputImageType > CastFilterTypeInverse;
    // auto filterInv = CastFilterTypeInverse::New();
    // // InputImageType i = *img->Clone();
    // filterInv->SetInput(mo);
    // filterInv->Update();
    // return filterInv->GetOutput();

}
Denoise::OtsuTypeP Denoise::BinaryOtsu(Denoise::RealTypeP img){
    typedef itk::OtsuThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> OtsuType;
    auto b = OtsuType::New();
    b->SetInput(img);
    b->SetInsideValue(1);
    b->SetOutsideValue(255);
//    b->Update();
    return b;
}
Denoise::OtsuTypeP Denoise::BinaryOtsu(Denoise::InputTypeP img){

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    filter->SetInput(img);
    filter->Update();
    return BinaryOtsu(filter->GetOutput());
}


Denoise::HuangTypeP Denoise::BinaryHuang(Denoise::RealTypeP img){
    typedef itk::HuangThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> HuangType;
    auto b = HuangType::New();
    b->SetInput(img);
    b->SetInsideValue(1);
    b->SetOutsideValue(255);
//    b->Update();
    return b;
}
Denoise::HuangTypeP Denoise::BinaryHuang(Denoise::InputTypeP img){

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    filter->SetInput(img);
    filter->Update();
    return BinaryHuang(filter->GetOutput());
}

Denoise::YenTypeP Denoise::BinaryYen(Denoise::RealTypeP img){
    typedef itk::YenThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> YenType;
    auto b = YenType::New();
    b->SetInput(img);
    b->SetInsideValue(1);
    b->SetOutsideValue(255);
//    b->Update();
    return b;
}
Denoise::YenTypeP Denoise::BinaryYen(Denoise::InputTypeP img){

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    filter->SetInput(img);
    filter->Update();
    return BinaryYen(filter->GetOutput());
}
Denoise::ShanbhagTypeP Denoise::BinaryShanbhag(Denoise::RealTypeP img){
    typedef itk::ShanbhagThresholdImageFilter<Denoise::RealImageType, Denoise::InputImageType> ShanbhagType;
    auto b = ShanbhagType::New();
    b->SetInput(img);
    b->SetInsideValue(1);
    b->SetOutsideValue(255);
//    b->Update();
    return b;
}
Denoise::ShanbhagTypeP Denoise::BinaryShanbhag(Denoise::InputTypeP img){

    typedef itk::CastImageFilter< Denoise::InputImageType, Denoise::RealImageType > CastFilterType;
    auto filter = CastFilterType::New();
    filter->SetInput(img);
    filter->Update();
    return BinaryShanbhag(filter->GetOutput());
}
