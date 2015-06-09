#ifndef itkFlatStructuringElementWithImageBridge_h
#define itkFlatStructuringElementWithImageBridge_h
#include "itkFlatStructuringElement.h"

namespace itk {

/** \class FlatStructuringElementWithImageBridge
 * \brief Extension to FlatStructuringElement to create elements from images, and to gather images from images. Based on 2006 dif, rejected because different compiler problems.
 *
 * For more detailed experiments with this mode, please refer to the results of the
 * test itkFlatStructuringElementTest.cxx or the wiki example.
 *
 * \ingroup ITKMathematicalMorphology
 *
 * \wiki
 * \wikiexample{Morphology/FlatStructuringElement,Erode a binary image using a flat (box) structuring element}
 * \endwiki
 * \wiki
 * \wikiexample{Morphology/FlatStructuringElementRadiusIsParametric,Generate structuring elements with accurate area}
 * \endwiki
 */

template<unsigned int VDimension, typename TImageType>
class FlatStructuringElementWithImageBridge:
    public FlatStructuringElement<VDimension> {
public:
  /** Standard class typedefs. */
  typedef FlatStructuringElementWithImageBridge< VDimension, TImageType > Self;
  typedef FlatStructuringElement< VDimension > Superclass;
  // typedef typename Superclass::Superclass     Neighborhood;

  /** Iterator typedef support. Note the naming is intentional, i.e.,
  * AllocatorType::iterator and AllocatorType::const_iterator, because the
  * allocator may be a vnl object or other type, which uses this form. */
  typedef typename Superclass::Iterator      Iterator;
  typedef typename Superclass::ConstIterator ConstIterator;

  /** Size and value typedef support. */
  typedef typename Superclass::SizeType      SizeType;
  typedef typename Superclass::OffsetType    OffsetType;

  /** Radius typedef support. */
  typedef typename Superclass::RadiusType RadiusType;

  /** External slice iterator type typedef support. */
  typedef typename Superclass::SliceIteratorType SliceIteratorType;

  /** External support for dimensionality. */
  itkStaticConstMacro(NeighborhoodDimension, unsigned int, VDimension);

  typedef Vector< float, VDimension > LType;
  typedef std::vector< LType >        DecompType;

  /** Typedefs related to the Image */
  typedef TImageType                     ImageType;
  typedef typename ImageType::Pointer    ImagePointer;
  typedef typename ImageType::PixelType  ImagePixelType;

  /** Default destructor. */
  virtual ~FlatStructuringElementWithImageBridge() {}

  /** Default constructor. */
  FlatStructuringElementWithImageBridge(): Superclass()
  {
  }
  /** Various constructors */
  static Self FromImage(
	                  const ImageType * image,
	                  ImagePixelType foreground =
	                      NumericTraits< ImagePixelType >::max() );

  // static Self FromImageUC(
  //         const UnsignedCharImageType * image,
  //         unsigned char foreground =
  //         NumericTraits< unsigned char >::max()
  //         );

  /** return an itk::Image from the structuring element. Background defaults to
   * NumericTraits< PixelType >::Zero and foreground to
   * NumericTraits< PixelType >::max()
   */
  ImagePointer GetImage(
          ImagePixelType foreground = NumericTraits< ImagePixelType >::max(),
          ImagePixelType background = NumericTraits< ImagePixelType >::Zero );

  // #<{(|* return an itk::Image< unsigned char, VDimension > from the structuring element
  //  * This method is there to be used from wrappers. From C++, you should prefer
  //  * the GetImage() method.
  //  |)}>#
  // UnsignedCharImagePointer GetImageUC( unsigned char foreground,
  //         unsigned char background );
  //
  // #<{(|* return an itk::Image< unsigned char, VDimension > from the structuring element
  //  * This method is there to be used from wrappers. From C++, you should prefer
  //  * the GetImage() method.
  //  |)}>#
  // UnsignedCharImagePointer GetImageUC();
};
}//itk namespace end
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFlatStructuringElementWithImageBridge.hxx"
#endif
#endif
