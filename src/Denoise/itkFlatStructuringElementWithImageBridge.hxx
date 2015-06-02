#ifndef itkFlatStructuringElementWithImageBridge_hxx
#define itkFlatStructuringElementWithImageBridge_hxx
#include "itkFlatStructuringElementWithImageBridge.h"
namespace itk{

template<unsigned int VDimension, class TImage>
FlatStructuringElementWithImageBridge<VDimension, TImage>
FlatStructuringElementWithImageBridge<VDimension, TImage>
::FromImage(const ImageType * image, ImagePixelType foreground)
{
    RadiusType size = image->GetLargestPossibleRegion().GetSize();
    Index< VDimension > centerIdx;

    for( unsigned int i=0; i<VDimension; i++ )
    {
        if( ( size[i] & 1 ) )
        {
            itk::ExceptionObject excp;
            excp.SetDescription("Size is not odd");
        }
        size[i] = size[i] / 2;
        centerIdx[i] = size[i];
    }
    FlatStructuringElementWithImageBridge res = FlatStructuringElementWithImageBridge<VDimension, TImage>();
    res.SetRadius( size );

    for( unsigned int j=0; j < res.Size(); j++ )
    {
        res[j] = image->GetPixel( centerIdx + res.GetOffset( j ) );
    }

    return res;
}


template<unsigned int VDimension, class TImage>
typename FlatStructuringElementWithImageBridge<VDimension, TImage>::ImagePointer
    FlatStructuringElementWithImageBridge<VDimension, TImage>
::GetImage(ImagePixelType foreground, ImagePixelType background)
{
    typename ImageType::Pointer image = ImageType::New();
    typename ImageType::RegionType region;
    RadiusType size = this->GetRadius();
    Index< VDimension > centerIdx;

    for( unsigned int i = 0; i < VDimension; i++ )
    {
        centerIdx[i] = size[i];
        size[i] = 2*size[i] + 1;
    }

    region.SetSize( size );
    image->SetRegions( region );
    image->Allocate();


    for(unsigned int j=0; j<this->Size(); j++ )
    {
        if( this->GetElement( j ) )
        {
            image->SetPixel( centerIdx+this->GetOffset( j ), foreground );
        }
        else
        {
            image->SetPixel( centerIdx+this->GetOffset( j ), background );
        }
    }

    return image;

}
}//ITK Namespace
#endif
