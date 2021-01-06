//
// With this program (using the ITK software library) you can upsample a given 3D volume of type float and apply deconvolution to increase detail.
//
//
// Parameters:
//              inFile
//              blurringKernel
//              outFile
//              sx,sy and sz : scaling parameters in three different directions
//              normalize kernel (yes/no)
//              number of iterations
//              Aplpha relaxation
//
// This software is is published under the GNU General Public License version 3(www.gnu.org/licenses/gpl-3.0.en.html)
// The author and University Medical Center Utrecht  make no representations about the
// suitability of this software for any purpose. It is provided "as is" without express or implied warranty.
//
// Rene Mandl version 1.0  2020
//

#include "itkFFTConvolutionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkLandweberDeconvolutionImageFilter.h"
#include "itkProjectedLandweberDeconvolutionImageFilter.h"
#include "itkRichardsonLucyDeconvolutionImageFilter.h"
#include "itkWienerDeconvolutionImageFilter.h"
#include "itkParametricBlindLeastSquaresDeconvolutionImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkMacro.h"
#include "itkMath.h"
#include <itkNiftiImageIO.h>
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include <metaCommand.h>

//
// getImageIO
//

itk::ImageIOBase::Pointer getImageIO(std::string input) {
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(input.c_str(), itk::ImageIOFactory::ReadMode);

  imageIO->SetFileName(input);
  imageIO->ReadImageInformation();

  return imageIO;
}


//
// MAIN
//

int main(int argc, char* argv[])
{
    //
    //  Parse program commands with Metacommand
    //
    
    MetaCommand command;
    
    command.AddField("infile", "name input volume file", MetaCommand::STRING,true);
    command.AddField("kernel", "name kernel volume file", MetaCommand::STRING,true);
    command.AddField("outfile", "name output volume file", MetaCommand::STRING,true);
    
    command.SetOption("noNormalize","n", false, "do not normalize the deconvolution kernel (default is normalize)");
    command.SetOptionLongTag("noNormalize","noNormalizeKernel");
    
    command.SetOption("subsampleFactor","s", false, "subsample factor (x y z) fo input volume");
    command.SetOptionLongTag("subsampleFactor","subsampleFactor");
    command.AddOptionField("subsampleFactor", "SX", MetaCommand::INT, false,"1");
    command.AddOptionField("subsampleFactor", "SY", MetaCommand::INT, false,"1");
    command.AddOptionField("subsampleFactor", "SZ", MetaCommand::INT, false,"1");

    command.SetOption("nrOfIterations","i", false, "number of iterations");
    command.SetOptionLongTag("nrOfIterations","nrOfIterations");
    command.AddOptionField("nrOfIterations", "Nr", MetaCommand::INT, true,"5");
    
    command.SetOption("Alpha","a", false, "relaxation factor alpha");
    command.SetOptionLongTag("Alpha","alpha");
    command.AddOptionField("Alpha", "alpha", MetaCommand::FLOAT, true,"0.05");
    
    command.SetOption("DeconvolutionType","t", false,"Deconvolution type: 1 = Landweber; 2 = projected Landweber; 3 = Richardson Lucy; 4 = Wiener; 5 = parametric blind least squares (FOR NOW ONLY '1' IS IMPLEMENTED) ");
    command.SetOptionLongTag("DeconvolutionType","deconvolutionType");
    command.AddOptionField("DeconvolutionType", "DeconvType", MetaCommand::INT, false,"1");
    
    if (!command.Parse(argc, argv)) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "With this program you can upsample a given 3D volume of type float and apply deconvolution to increase detail." << std::endl;
        std::cerr << std::endl;
        std::cerr << "USAGE:" << std::endl;
        std::cerr << "itkLandweberDeconvolution <inputVolume>.nii <kernel>.nii <outputVolume>.nii" << std::endl;
        return EXIT_FAILURE;
    }

    std::string inFileName = command.GetValueAsString("infile");
    std::string kernelFileName = command.GetValueAsString("kernel");
    std::string outFileName = command.GetValueAsString("outfile");
    
    bool noNormalizeKernel = command.GetValueAsBool("noNormalize");
    int subsampleFactorX   = command.GetValueAsInt("subsampleFactor","SX");
    int subsampleFactorY   = command.GetValueAsInt("subsampleFactor","SY");
    int subsampleFactorZ   = command.GetValueAsInt("subsampleFactor","SZ");
    int nrOfIterations     = command.GetValueAsInt("nrOfIterations","Nr");
    float alpha            = command.GetValueAsFloat("Alpha","alpha");
    int deconvolutionType  = command.GetValueAsInt("DeconvolutionType","DeconvType");
    
    
    // EOF command parsing
    
    const unsigned int Dimension = 3;
    
    using inputVoxelType = float; // for now only operate on float images
    using ImageType = itk::Image<inputVoxelType, Dimension>;
    using WriterType = itk::ImageFileWriter<ImageType>;
    using LandWeberDeconvolutionFilterType = itk::LandweberDeconvolutionImageFilter< ImageType >;
    using DeconvolutionFilterType = itk::LandweberDeconvolutionImageFilter< ImageType >;
    
    itk::ImageIOBase::Pointer volumeIO = getImageIO(inFileName);
    itk::ImageIOBase::Pointer filterIO = getImageIO(kernelFileName);

    using inVolReader = itk::ImageFileReader<ImageType>;

    inVolReader::Pointer in_volReader = inVolReader::New();
    in_volReader->SetFileName(volumeIO->GetFileName());
    in_volReader->Update();

    using inFilterReader = itk::ImageFileReader<ImageType>;

    // resampler
    const ImageType::SpacingType& inputSpacing = in_volReader->GetOutput()->GetSpacing();
    ImageType::SpacingType outputSpacing;
    outputSpacing[0] = inputSpacing[0] / subsampleFactorX;
    outputSpacing[1] = inputSpacing[1] / subsampleFactorY;
    outputSpacing[2] = inputSpacing[2] / subsampleFactorZ;
    
    ImageType::SizeType inputSize=in_volReader->GetOutput()->GetLargestPossibleRegion().GetSize();
    ImageType::SizeType outputSize;
    outputSize[0] = inputSize[0] * subsampleFactorX;
    outputSize[1] = inputSize[1] * subsampleFactorY;
    outputSize[2] = inputSize[2] * subsampleFactorZ;
    
    using filterType=itk::ResampleImageFilter<ImageType, ImageType>;
    filterType::Pointer resampler = filterType::New();
    using transformType=itk::IdentityTransform<double, Dimension>;
    transformType::Pointer transform = transformType::New();
    transform->SetIdentity();
    resampler->SetTransform(transform);
    resampler->SetOutputSpacing(outputSpacing);
    resampler->SetDefaultPixelValue(0); // should never occur
    resampler->SetOutputDirection(in_volReader->GetOutput()->GetDirection());
    resampler->SetOutputOrigin(in_volReader->GetOutput()->GetOrigin());
    resampler->SetSize(outputSize);
    
    resampler->SetInput(in_volReader->GetOutput());
    
    // the deblurring kernel
    inFilterReader::Pointer in_filterReader = inFilterReader::New();
    in_filterReader->SetFileName(filterIO->GetFileName());
    in_filterReader->Update();

    // To implement if necessary: different types of the deconvolution algorithm
    // if (1 == deconvolutionType) {DeconvolutionFilterType::Pointer deconvolutionFilter = LandWeberDeconvolutionFilterType::New();}
    // if (2 == deconvolutionType) {DeconvolutionFilterType::Pointer deconvolutionFilter ProjectedLandWeberDeconvolutionFilterType::New();}
    // if (3 == deconvolutionType) {DeconvolutionFilterType::Pointer deconvolutionFilter = RichardsonLucyDeconvolutionFilterType::New();}
    // if (4 == deconvolutionType) {DeconvolutionFilterType::Pointer deconvolutionFilter = WienerDeconvolutionFilterType::New();}
    // if (5 == deconvolutionType) {DeconvolutionFilterType::Pointer deconvolutionFilter = ParametricBlindLeastSquaresDeconvolutionFilterType::New();}

    // the statement below will be replaced by the option for the various different filter types
    DeconvolutionFilterType::Pointer deconvolutionFilter = DeconvolutionFilterType::New();
    
    deconvolutionFilter->SetInput(resampler->GetOutput());
    deconvolutionFilter->SetKernelImage(in_filterReader->GetOutput());
    
    deconvolutionFilter->NormalizeOn();
    if (noNormalizeKernel)deconvolutionFilter->NormalizeOff();
    
    deconvolutionFilter->SetAlpha(alpha);
    
    if (itk::Math::NotExactlyEquals(deconvolutionFilter->GetAlpha(), alpha)) {
        std::cerr << "Set/GetAlpha() test failed." << std::endl;
        return EXIT_FAILURE;
    }

    deconvolutionFilter->SetNumberOfIterations(nrOfIterations);

    // Write out the deconvolution result
    try {
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName(outFileName);
        writer->SetInput(deconvolutionFilter->GetOutput());
        writer->Update();
    } catch ( itk::ExceptionObject & e ) {
        std::cerr << "Unexpected exception caught when writing deconvolution image: " << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
