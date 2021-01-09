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
#include "itkTikhonovDeconvolutionImageFilter.h"
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
    
    command.SetOption("subsampleFactor","s", false, "subsample factor (x y z) fo input volume (should be integer and 1 or larger)");
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
    
    command.SetOption("DeconvolutionType","t", false,"Deconvolution type: 1 = Landweber; 2 = projected Landweber; 3 = Richardson Lucy; 4 = Wiener; 5 = Tikhonov");
    command.SetOptionLongTag("DeconvolutionType","deconvolutionType");
    command.AddOptionField("DeconvolutionType", "DeconvType", MetaCommand::INT, false,"1");
    
    if (!command.Parse(argc, argv)) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "With this program you can upsample a given 3D volume of type float and apply deconvolution to increase detail." << std::endl;
        std::cerr << std::endl;
        std::cerr << "USAGE:" << std::endl;
        std::cerr << "itkDeconvolution <inputVolume>.nii <kernel>.nii <outputVolume>.nii" << std::endl;
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
    
    if (subsampleFactorX < 1 || subsampleFactorY < 1 || subsampleFactorZ < 1) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "ERROR: subsample factors should be integer and equal or larger than 1" << std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }
    
    if ( deconvolutionType < 1 || deconvolutionType > 5 ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "ERROR: deconvolution type should be 1-5" << std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }
    
    if ( nrOfIterations < 0 ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "ERROR: nr of iterations should be equal or larger than zero" << std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }
    
    // EOF command parsing
    
    const unsigned int Dimension = 3;
    
    using inputVoxelTypeFloat = float;
    using inputVoxelTypeInt = int;
    
    using ImageType = itk::Image<inputVoxelTypeFloat, Dimension>;
    using WriterType = itk::ImageFileWriter<ImageType>;
    
    using LandWeberDeconvolutionFilterType = itk::LandweberDeconvolutionImageFilter< ImageType >;
    using ProjectedLandweberDeconvolutionFilterType = itk::ProjectedLandweberDeconvolutionImageFilter< ImageType >;
    using RichardsonLucyDeconvolutionFilterType = itk::RichardsonLucyDeconvolutionImageFilter< ImageType >;
    using WienerDeconvolutionFilterType = itk::WienerDeconvolutionImageFilter< ImageType >;
    using TikhonovDeconvolutionImageFilterType = itk::TikhonovDeconvolutionImageFilter< ImageType >;
    
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

    
    LandWeberDeconvolutionFilterType::Pointer LandWeberDeconvolutionFilter = LandWeberDeconvolutionFilterType::New();
    ProjectedLandweberDeconvolutionFilterType::Pointer ProjectedLandWeberDeconvolutionFilter = ProjectedLandweberDeconvolutionFilterType::New();
    RichardsonLucyDeconvolutionFilterType::Pointer RichardsonLucyDeconvolutionFilter = RichardsonLucyDeconvolutionFilterType::New();
    WienerDeconvolutionFilterType::Pointer WienerDeconvolutionFilter = WienerDeconvolutionFilterType::New();
    TikhonovDeconvolutionImageFilterType::Pointer TikhonovDeconvolutionFilter = TikhonovDeconvolutionImageFilterType::New();
    
    switch (deconvolutionType) {
        case 1:
            LandWeberDeconvolutionFilter->SetInput(resampler->GetOutput());
            LandWeberDeconvolutionFilter->SetKernelImage(in_filterReader->GetOutput());
            
            LandWeberDeconvolutionFilter->NormalizeOn();
            if (noNormalizeKernel)LandWeberDeconvolutionFilter->NormalizeOff();
            
            LandWeberDeconvolutionFilter->SetAlpha(alpha);
            
            if (itk::Math::NotExactlyEquals(LandWeberDeconvolutionFilter->GetAlpha(), alpha)) {
                std::cerr << "Set/GetAlpha() test failed." << std::endl;
                return EXIT_FAILURE;
            }

            LandWeberDeconvolutionFilter->SetNumberOfIterations(nrOfIterations);
            LandWeberDeconvolutionFilter->SetOutputRegionModeToSame();
            LandWeberDeconvolutionFilter->Update();
            
            // Write out the deconvolution result
            try {
                WriterType::Pointer writer = WriterType::New();
                writer->SetFileName(outFileName);
                writer->SetInput(LandWeberDeconvolutionFilter->GetOutput());
                writer->Update();
            } catch ( itk::ExceptionObject & e ) {
                std::cerr << "Unexpected exception caught when writing deconvoluted image: " << std::endl;
                return EXIT_FAILURE;
            }
            break;
            
        case 2:
            ProjectedLandWeberDeconvolutionFilter->SetInput(resampler->GetOutput());
            ProjectedLandWeberDeconvolutionFilter->SetKernelImage(in_filterReader->GetOutput());
            
            ProjectedLandWeberDeconvolutionFilter->NormalizeOn();
            if (noNormalizeKernel)ProjectedLandWeberDeconvolutionFilter->NormalizeOff();
            
            ProjectedLandWeberDeconvolutionFilter->SetAlpha(alpha);
            
            if (itk::Math::NotExactlyEquals(ProjectedLandWeberDeconvolutionFilter->GetAlpha(), alpha)) {
                std::cerr << "Set/GetAlpha() test failed." << std::endl;
                return EXIT_FAILURE;
            }

            ProjectedLandWeberDeconvolutionFilter->SetNumberOfIterations(nrOfIterations);
            ProjectedLandWeberDeconvolutionFilter->SetOutputRegionModeToSame();
            ProjectedLandWeberDeconvolutionFilter->Update();
            
            // Write out the deconvolution result
            try {
                WriterType::Pointer writer = WriterType::New();
                writer->SetFileName(outFileName);
                writer->SetInput(ProjectedLandWeberDeconvolutionFilter->GetOutput());
                writer->Update();
            } catch ( itk::ExceptionObject & e ) {
                std::cerr << "Unexpected exception caught when writing deconvoluted image: " << std::endl;
                return EXIT_FAILURE;
            }
            break;
        case 3:
            RichardsonLucyDeconvolutionFilter->SetInput(resampler->GetOutput());
            RichardsonLucyDeconvolutionFilter->SetKernelImage(in_filterReader->GetOutput());
            
            RichardsonLucyDeconvolutionFilter->NormalizeOn();
            if (noNormalizeKernel)RichardsonLucyDeconvolutionFilter->NormalizeOff();

            RichardsonLucyDeconvolutionFilter->SetNumberOfIterations(nrOfIterations);
            RichardsonLucyDeconvolutionFilter->SetOutputRegionModeToSame();
            RichardsonLucyDeconvolutionFilter->Update();
            
            // Write out the deconvolution result
            try {
                WriterType::Pointer writer = WriterType::New();
                writer->SetFileName(outFileName);
                writer->SetInput(RichardsonLucyDeconvolutionFilter->GetOutput());
                writer->Update();
            } catch ( itk::ExceptionObject & e ) {
                std::cerr << "Unexpected exception caught when writing deconvoluted image: " << std::endl;
                return EXIT_FAILURE;
            }
            break;
        case 4:
            WienerDeconvolutionFilter->SetInput(resampler->GetOutput());
            WienerDeconvolutionFilter->SetKernelImage(in_filterReader->GetOutput());
            
            WienerDeconvolutionFilter->NormalizeOn();
            if (noNormalizeKernel)WienerDeconvolutionFilter->NormalizeOff();
            
            WienerDeconvolutionFilter->SetOutputRegionModeToSame();
            WienerDeconvolutionFilter->Update();
            
            // Write out the deconvolution result
            try {
                WriterType::Pointer writer = WriterType::New();
                writer->SetFileName(outFileName);
                writer->SetInput(WienerDeconvolutionFilter->GetOutput());
                writer->Update();
            } catch ( itk::ExceptionObject & e ) {
                std::cerr << "Unexpected exception caught when writing deconvoluted image: " << std::endl;
                return EXIT_FAILURE;
            }
            break;
            
        case 5:
            TikhonovDeconvolutionFilter->SetInput(resampler->GetOutput());
            TikhonovDeconvolutionFilter->SetKernelImage(in_filterReader->GetOutput());
            
            TikhonovDeconvolutionFilter->NormalizeOn();
            if (noNormalizeKernel) TikhonovDeconvolutionFilter->NormalizeOff();
            
            TikhonovDeconvolutionFilter->SetOutputRegionModeToSame();
            TikhonovDeconvolutionFilter->Update();
            
            // Write out the deconvolution result
            try {
                WriterType::Pointer writer = WriterType::New();
                writer->SetFileName(outFileName);
                writer->SetInput(TikhonovDeconvolutionFilter->GetOutput());
                writer->Update();
            } catch ( itk::ExceptionObject & e ) {
                std::cerr << "Unexpected exception caught when writing deconvoluted image: " << std::endl;
                return EXIT_FAILURE;
            }
            break;
        default:
            std::cerr << std::endl << std::endl;
            std::cerr << "Unknown dconvolution type!" << std::endl;
            return EXIT_FAILURE;
            break;
    }

    
    return EXIT_SUCCESS;
}
