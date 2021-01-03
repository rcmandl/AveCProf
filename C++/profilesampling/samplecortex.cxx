//
// samplecortex.cxx
//
// This program samples cortex values at different depths from a nifti file.
//
//
// USAGE: <pialCoordinatesFile> <wmCoordinatesFile> <dataToSample.mnc>
//
//        Coordinate files (ASCII) should contain an x y z coordinate at each line. It is assumed that there is a point-to-point correspondence between
//        the coordinate files
//        Each line in output starts with segmentnr and point coordinates of surfaces, then the measurement values and ends with the relative resolution value.
//
//
// Note:
// There is something strange when reading in the nifti files. For some reason the x and y coordinates from scanner space need to be negated when calling inVolume->TransformPhysicalPointToIndex()
// WHen comparing the directionality matrix and origin with the matrices from say fslhd, you can see that indeed the x and y are negated. (it appears not to be related to RAS vs LPS)
// For now an option is added allowing to mirror x and y values but this need to be further clarified.
//
// This software is is published under the GNU General Public License version 3(www.gnu.org/licenses/gpl-3.0.en.html)
// The author and University Medical Center Utrecht  make no representations about the
// suitability of this software for any purpose. It is provided "as is" without express or implied warranty.
//
// Rene Mandl version 1.0  2020

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMacro.h"
#include "itkMath.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkOrientImageFilter.h"
#include <itkNiftiImageIO.h>
#include <metaCommand.h>
#include <iostream>
#include <cstdio>

const unsigned int Dimension = 3; // only operate on 3D volumes

using inputVoxelType = float; // for now only operate on float images
using ImageType = itk::Image<inputVoxelType, Dimension>;
using linearInterpolatorType = itk::LinearInterpolateImageFunction<ImageType, inputVoxelType>;// Needed if we want to do linear interpolation instead of nearest neighbour
using bSplineInterpolatorType = itk::BSplineInterpolateImageFunction<ImageType, inputVoxelType>;// Needed if we want to do bSpline interpolation instead of nearest neighbour
using orientImageFilterType = itk::OrientImageFilter<ImageType, ImageType>; // to change the coordidnate system for the nifti files
//
// getImageIO
//

itk::ImageIOBase::Pointer getImageIO(std::string input) {
    itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(input.c_str(), itk::ImageIOFactory::ReadMode);
    
    if ( NULL == imageIO ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "*** ERROR: Could not read file: " << input << std::endl;
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }
    imageIO->SetFileName(input);
    imageIO->ReadImageInformation();
    
    return imageIO;
}


//
// sampleAndWriteLineWithFlagsNN
//
// Perform sampling using nearest neighbour and add to the flagVolume

void sampleAndWriteLineWithFlagsNN( ImageType::Pointer inVolume, ImageType::Pointer flagVolume, float pialCoordinateX, float pialCoordinateY, float pialCoordinateZ, float wmCoordinateX, float wmCoordinateY, float wmCoordinateZ, unsigned int nrOfSteps, float extentFactor, long lineNr, float SCL ) {
    
    unsigned int i;
    using PointType = itk::Point<double, Dimension >;
    PointType point;
    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;
    
    float length, extPartLength, normX, normY, normZ, xStep, yStep, zStep;
    float extPialCoordinateX, extPialCoordinateY, extPialCoordinateZ, extWmCoordinateX, extWmCoordinateY, extWmCoordinateZ;
    
    // length of the original segment
    length = sqrt( ( pialCoordinateX - wmCoordinateX ) * ( pialCoordinateX - wmCoordinateX ) + ( pialCoordinateY - wmCoordinateY ) * ( pialCoordinateY - wmCoordinateY ) +( pialCoordinateZ - wmCoordinateZ ) * ( pialCoordinateZ - wmCoordinateZ ) );
    
    if ( 0 >= length ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "*** WARNING: profile length <= zero!!" << std::endl;
        std::cerr << std::endl;
    } else {
        // length of the extention at both sides
        extPartLength = (( length * extentFactor ) - length) / 2;
        
        normX = ( pialCoordinateX - wmCoordinateX ) / length;
        normY = ( pialCoordinateY - wmCoordinateY ) / length;
        normZ = ( pialCoordinateZ - wmCoordinateZ ) / length;
        
        // compute new pial point including extention for pial boundary
        extPialCoordinateX = pialCoordinateX + (normX * extPartLength);
        extPialCoordinateY = pialCoordinateY + (normY * extPartLength);
        extPialCoordinateZ = pialCoordinateZ + (normZ * extPartLength);
        
        // for extention for wm boundary multiply normX etc by -1
        extWmCoordinateX = wmCoordinateX + (-1 * normX * extPartLength);
        extWmCoordinateY = wmCoordinateY + (-1 * normY * extPartLength);
        extWmCoordinateZ = wmCoordinateZ + (-1 * normZ * extPartLength);
        
        // step computed from wm towards pial surface
        // note that we want a certain nr of steps (not segments) so we go from 0 to (nrOfSteps - 1)
        xStep = ( extPialCoordinateX - extWmCoordinateX ) / ( nrOfSteps - 1);
        yStep = ( extPialCoordinateY - extWmCoordinateY ) / ( nrOfSteps - 1);
        zStep = ( extPialCoordinateZ - extWmCoordinateZ ) / ( nrOfSteps - 1);
        
        point[0]=extWmCoordinateX;
        point[1]=extWmCoordinateY;
        point[2]=extWmCoordinateZ;
        
        
        std::cout << lineNr << " " << extPialCoordinateX << " " << extPialCoordinateY << " " << extPialCoordinateZ << " " << extWmCoordinateX << " " << extWmCoordinateY << " " << extWmCoordinateZ << " ";
        
        for (i = 0; i < nrOfSteps; i++) {
            
            const bool isInside = inVolume->TransformPhysicalPointToIndex(point,pixelIndex); // test if coordinate falls within image
            
            if ( isInside ) {
                pixelValue = inVolume->GetPixel(pixelIndex);
                std::cout << pixelValue << " ";
                //
                // Now increment the value of the corresponding voxel in the flagVolume
                //
                pixelValue= flagVolume->GetPixel(pixelIndex);
                pixelValue++;
                flagVolume->SetPixel(pixelIndex,pixelValue);
            } else {
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "*** ERROR: coordinate of sample point is outside volume; should not occur! (perhaps -negateXY option is needed?)" << std::endl;
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
            }
            
            point[0]+=xStep;
            point[1]+=yStep;
            point[2]+=zStep;
        }
        std::cout << length << " " << SCL << std::endl;
    }
}


//
// sampleAndWriteLineNN
//
// Perform sampling using nearest neighbour

void sampleAndWriteLineNN( ImageType::Pointer inVolume, float pialCoordinateX, float pialCoordinateY, float pialCoordinateZ, float wmCoordinateX, float wmCoordinateY, float wmCoordinateZ, unsigned int nrOfSteps, float extentFactor, long lineNr, float SCL ) {
    
    unsigned int i;
    using PointType = itk::Point<double, Dimension >;
    PointType point;
    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;
    
    float length, extPartLength, normX, normY, normZ, xStep, yStep, zStep;
    float extPialCoordinateX, extPialCoordinateY, extPialCoordinateZ, extWmCoordinateX, extWmCoordinateY, extWmCoordinateZ;
    
    // length of the original segment
    length = sqrt( ( pialCoordinateX - wmCoordinateX ) * ( pialCoordinateX - wmCoordinateX ) + ( pialCoordinateY - wmCoordinateY ) * ( pialCoordinateY - wmCoordinateY ) +( pialCoordinateZ - wmCoordinateZ ) * ( pialCoordinateZ - wmCoordinateZ ) );
    
    if ( 0 >= length ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "*** WARNING: profile length <= zero!!" << std::endl;
        std::cerr << std::endl;
    } else {
        // length of the extention at both sides
        extPartLength = (( length * extentFactor ) - length) / 2;
    
        normX = ( pialCoordinateX - wmCoordinateX ) / length;
        normY = ( pialCoordinateY - wmCoordinateY ) / length;
        normZ = ( pialCoordinateZ - wmCoordinateZ ) / length;
    
        // compute new pial point including extention for pial boundary
        extPialCoordinateX = pialCoordinateX + (normX * extPartLength);
        extPialCoordinateY = pialCoordinateY + (normY * extPartLength);
        extPialCoordinateZ = pialCoordinateZ + (normZ * extPartLength);
    
        // for extention for wm boundary multiply normX etc by -1
        extWmCoordinateX = wmCoordinateX + (-1 * normX * extPartLength);
        extWmCoordinateY = wmCoordinateY + (-1 * normY * extPartLength);
        extWmCoordinateZ = wmCoordinateZ + (-1 * normZ * extPartLength);
    
        // step computed from wm towards pial surface
        // note that we want a certain nr of steps (not segments) so we go from 0 to (nrOfSteps - 1)
        xStep = ( extPialCoordinateX - extWmCoordinateX ) / ( nrOfSteps - 1);
        yStep = ( extPialCoordinateY - extWmCoordinateY ) / ( nrOfSteps - 1);
        zStep = ( extPialCoordinateZ - extWmCoordinateZ ) / ( nrOfSteps - 1);
    
        point[0]=extWmCoordinateX;
        point[1]=extWmCoordinateY;
        point[2]=extWmCoordinateZ;
    
        std::cout << lineNr << " " << extPialCoordinateX << " " << extPialCoordinateY << " " << extPialCoordinateZ << " " << extWmCoordinateX << " " << extWmCoordinateY << " " << extWmCoordinateZ << " ";
        
        for (i = 0; i < nrOfSteps; i++) {
            
            const bool isInside = inVolume->TransformPhysicalPointToIndex(point,pixelIndex); // test if coordinate falls within image
            
            if ( isInside ) {
                pixelValue = inVolume->GetPixel(pixelIndex);
                std::cout << pixelValue << " ";
            } else {
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "*** ERROR: coordinate of sample point is outside volume (should not occur!)" << std::endl;
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
            }

            point[0]+=xStep;
            point[1]+=yStep;
            point[2]+=zStep;
        }
        std::cout << length << " " << SCL << std::endl;
    }
}



//
// sampleAndWriteLineWithFlagsLI
//
// Perform sampling using linear interpolation and add to flagVolume

void sampleAndWriteLineWithFlagsLI( ImageType::Pointer inVolume, linearInterpolatorType::Pointer inter, ImageType::Pointer flagVolume, float pialCoordinateX, float pialCoordinateY, float pialCoordinateZ, float wmCoordinateX, float wmCoordinateY, float wmCoordinateZ, unsigned int nrOfSteps, float extentFactor, long lineNr, float SCL ) {
    
    unsigned int i;
    using PointType = itk::Point<double, Dimension >;
    PointType point;
    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;
    
    float length, extPartLength, normX, normY, normZ, xStep, yStep, zStep;
    float extPialCoordinateX, extPialCoordinateY, extPialCoordinateZ, extWmCoordinateX, extWmCoordinateY, extWmCoordinateZ;
    
    // length of the original segment
    length = sqrt( ( pialCoordinateX - wmCoordinateX ) * ( pialCoordinateX - wmCoordinateX ) + ( pialCoordinateY - wmCoordinateY ) * ( pialCoordinateY - wmCoordinateY ) +( pialCoordinateZ - wmCoordinateZ ) * ( pialCoordinateZ - wmCoordinateZ ) );
    
    if ( 0 >= length ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "*** WARNING: profile length <= zero!!" << std::endl;
        std::cerr << std::endl;
    } else {
        // length of the extention at both sides
        extPartLength = (( length * extentFactor ) - length) / 2;
        
        normX = ( pialCoordinateX - wmCoordinateX ) / length;
        normY = ( pialCoordinateY - wmCoordinateY ) / length;
        normZ = ( pialCoordinateZ - wmCoordinateZ ) / length;
        
        // compute new pial point including extention for pial boundary
        extPialCoordinateX = pialCoordinateX + (normX * extPartLength);
        extPialCoordinateY = pialCoordinateY + (normY * extPartLength);
        extPialCoordinateZ = pialCoordinateZ + (normZ * extPartLength);
        
        // for extention for wm boundary multiply normX etc by -1
        extWmCoordinateX = wmCoordinateX + (-1 * normX * extPartLength);
        extWmCoordinateY = wmCoordinateY + (-1 * normY * extPartLength);
        extWmCoordinateZ = wmCoordinateZ + (-1 * normZ * extPartLength);
        
        // step computed from wm towards pial surface
        // note that we want a certain nr of steps (not segments) so we go from 0 to (nrOfSteps - 1)
        xStep = ( extPialCoordinateX - extWmCoordinateX ) / ( nrOfSteps - 1);
        yStep = ( extPialCoordinateY - extWmCoordinateY ) / ( nrOfSteps - 1);
        zStep = ( extPialCoordinateZ - extWmCoordinateZ ) / ( nrOfSteps - 1);
        
        point[0]=extWmCoordinateX;
        point[1]=extWmCoordinateY;
        point[2]=extWmCoordinateZ;
        
        std::cout << lineNr << " " << extPialCoordinateX << " " << extPialCoordinateY << " " << extPialCoordinateZ << " " << extWmCoordinateX << " " << extWmCoordinateY << " " << extWmCoordinateZ << " ";
        
        for (i = 0; i < nrOfSteps; i++) {
            
            const bool isInside = inVolume->TransformPhysicalPointToIndex(point,pixelIndex); // test if coordinate falls within image

            if ( isInside ) {
                pixelValue = inter->Evaluate(point); // Now use linear interpolation:
                std::cout << pixelValue << " ";
                //
                // Now increment the value of the corresponding voxel in the flagVolume
                //
                pixelValue= flagVolume->GetPixel(pixelIndex);
                pixelValue++;
                flagVolume->SetPixel(pixelIndex,pixelValue);
            } else {
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "*** ERROR: coordinate of sample point is outside volume (should not occur!)" << std::endl;
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
            }
            
            point[0]+=xStep;
            point[1]+=yStep;
            point[2]+=zStep;
        }
        std::cout << length << " " << SCL << std::endl;
    }
}




//
// sampleAndWriteLineLI
//
// Perform sampling using linear interpolation
//

void sampleAndWriteLineLI( ImageType::Pointer inVolume, linearInterpolatorType::Pointer inter, float pialCoordinateX, float pialCoordinateY, float pialCoordinateZ, float wmCoordinateX, float wmCoordinateY, float wmCoordinateZ, unsigned int nrOfSteps, float extentFactor, long lineNr, float SCL  ) {
    
    unsigned int i;
    using PointType = itk::Point<double, Dimension >;
    PointType point;
    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;
    
    float length, extPartLength, normX, normY, normZ, xStep, yStep, zStep;
    float extPialCoordinateX, extPialCoordinateY, extPialCoordinateZ, extWmCoordinateX, extWmCoordinateY, extWmCoordinateZ;
    
    // length of the original segment
    length = sqrt( ( pialCoordinateX - wmCoordinateX ) * ( pialCoordinateX - wmCoordinateX ) + ( pialCoordinateY - wmCoordinateY ) * ( pialCoordinateY - wmCoordinateY ) +( pialCoordinateZ - wmCoordinateZ ) * ( pialCoordinateZ - wmCoordinateZ ) );
    
    if ( 0 >= length ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "*** WARNING: profile length <= zero!!" << std::endl;
        std::cerr << std::endl;
    } else {
        // length of the extention at both sides
        extPartLength = (( length * extentFactor ) - length) / 2;
        
        normX = ( pialCoordinateX - wmCoordinateX ) / length;
        normY = ( pialCoordinateY - wmCoordinateY ) / length;
        normZ = ( pialCoordinateZ - wmCoordinateZ ) / length;
        
        // compute new pial point including extention for pial boundary
        extPialCoordinateX = pialCoordinateX + (normX * extPartLength);
        extPialCoordinateY = pialCoordinateY + (normY * extPartLength);
        extPialCoordinateZ = pialCoordinateZ + (normZ * extPartLength);
        
        // for extention for wm boundary multiply normX etc by -1
        extWmCoordinateX = wmCoordinateX + (-1 * normX * extPartLength);
        extWmCoordinateY = wmCoordinateY + (-1 * normY * extPartLength);
        extWmCoordinateZ = wmCoordinateZ + (-1 * normZ * extPartLength);
        
        // step computed from wm towards pial surface
        // note that we want a certain nr of steps (not segments) so we go from 0 to (nrOfSteps - 1)
        xStep = ( extPialCoordinateX - extWmCoordinateX ) / ( nrOfSteps - 1);
        yStep = ( extPialCoordinateY - extWmCoordinateY ) / ( nrOfSteps - 1);
        zStep = ( extPialCoordinateZ - extWmCoordinateZ ) / ( nrOfSteps - 1);
        
        point[0]=extWmCoordinateX;
        point[1]=extWmCoordinateY;
        point[2]=extWmCoordinateZ;
        
        
        std::cout << lineNr << " " << extPialCoordinateX << " " << extPialCoordinateY << " " << extPialCoordinateZ << " " << extWmCoordinateX << " " << extWmCoordinateY << " " << extWmCoordinateZ << " ";
        
        
        for (i = 0; i < nrOfSteps; i++) {
            
            const bool isInside = inVolume->TransformPhysicalPointToIndex(point, pixelIndex); // test if coordinate falls within image
            
            if ( isInside ) {
                pixelValue = inter->Evaluate(point); // Now use linear interpolation
                std::cout << pixelValue << " ";
            } else {
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "*** ERROR: coordinate of sample point is outside volume (should not occur!)" << std::endl;
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
            }
            
            point[0]+=xStep;
            point[1]+=yStep;
            point[2]+=zStep;
        }
        std::cout << length << " " << SCL << std::endl;
    }
}


//
// sampleAndWriteLineWithFlagsBS
//
// Perform sampling using linear interpolation and add to flagVolume

void sampleAndWriteLineWithFlagsBS( ImageType::Pointer inVolume, bSplineInterpolatorType::Pointer inter, ImageType::Pointer flagVolume, float pialCoordinateX, float pialCoordinateY, float pialCoordinateZ, float wmCoordinateX, float wmCoordinateY, float wmCoordinateZ, unsigned int nrOfSteps, float extentFactor, long lineNr, float SCL ) {
    
    unsigned int i;
    using PointType = itk::Point<double, Dimension >;
    PointType point;
    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;
    
    float length, extPartLength, normX, normY, normZ, xStep, yStep, zStep;
    float extPialCoordinateX, extPialCoordinateY, extPialCoordinateZ, extWmCoordinateX, extWmCoordinateY, extWmCoordinateZ;
    
    // length of the original segment
    length = sqrt( ( pialCoordinateX - wmCoordinateX ) * ( pialCoordinateX - wmCoordinateX ) + ( pialCoordinateY - wmCoordinateY ) * ( pialCoordinateY - wmCoordinateY ) +( pialCoordinateZ - wmCoordinateZ ) * ( pialCoordinateZ - wmCoordinateZ ) );
    
    if ( 0 >= length ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "*** WARNING: profile length <= zero!!" << std::endl;
        std::cerr << std::endl;
    } else {
        // length of the extention at both sides
        extPartLength = (( length * extentFactor ) - length) / 2;
        
        normX = ( pialCoordinateX - wmCoordinateX ) / length;
        normY = ( pialCoordinateY - wmCoordinateY ) / length;
        normZ = ( pialCoordinateZ - wmCoordinateZ ) / length;
        
        // compute new pial point including extention for pial boundary
        extPialCoordinateX = pialCoordinateX + (normX * extPartLength);
        extPialCoordinateY = pialCoordinateY + (normY * extPartLength);
        extPialCoordinateZ = pialCoordinateZ + (normZ * extPartLength);
        
        // for extention for wm boundary multiply normX etc by -1
        extWmCoordinateX = wmCoordinateX + (-1 * normX * extPartLength);
        extWmCoordinateY = wmCoordinateY + (-1 * normY * extPartLength);
        extWmCoordinateZ = wmCoordinateZ + (-1 * normZ * extPartLength);
        
        // step computed from wm towards pial surface
        // note that we want a certain nr of steps (not segments) so we go from 0 to (nrOfSteps - 1)
        xStep = ( extPialCoordinateX - extWmCoordinateX ) / ( nrOfSteps - 1);
        yStep = ( extPialCoordinateY - extWmCoordinateY ) / ( nrOfSteps - 1);
        zStep = ( extPialCoordinateZ - extWmCoordinateZ ) / ( nrOfSteps - 1);
        
        point[0]=extWmCoordinateX;
        point[1]=extWmCoordinateY;
        point[2]=extWmCoordinateZ;
        
        std::cout << lineNr << " " << extPialCoordinateX << " " << extPialCoordinateY << " " << extPialCoordinateZ << " " << extWmCoordinateX << " " << extWmCoordinateY << " " << extWmCoordinateZ << " ";
        
        for (i = 0; i < nrOfSteps; i++) {
            
            const bool isInside = inVolume->TransformPhysicalPointToIndex(point,pixelIndex); // test if coordinate falls within image
            
            if ( isInside ) {
                pixelValue = inter->Evaluate(point); // Now use linear interpolation:
                std::cout << pixelValue << " ";
                //
                // Now increment the value of the corresponding voxel in the flagVolume
                //
                pixelValue= flagVolume->GetPixel(pixelIndex);
                pixelValue++;
                flagVolume->SetPixel(pixelIndex,pixelValue);
            } else {
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "*** ERROR: coordinate of sample point is outside volume (should not occur!)" << std::endl;
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
            }
            
            point[0]+=xStep;
            point[1]+=yStep;
            point[2]+=zStep;
        }
        std::cout << length << " " << SCL << std::endl;
    }
}




//
// sampleAndWriteLineBS
//
// Perform sampling using linear interpolation
//

void sampleAndWriteLineBS( ImageType::Pointer inVolume, bSplineInterpolatorType::Pointer inter, float pialCoordinateX, float pialCoordinateY, float pialCoordinateZ, float wmCoordinateX, float wmCoordinateY, float wmCoordinateZ, unsigned int nrOfSteps, float extentFactor, long lineNr, float SCL ) {
    
    unsigned int i;
    using PointType = itk::Point<double, Dimension >;
    PointType point;
    ImageType::IndexType pixelIndex;
    ImageType::PixelType pixelValue;
    
    float length, extPartLength, normX, normY, normZ, xStep, yStep, zStep;
    float extPialCoordinateX, extPialCoordinateY, extPialCoordinateZ, extWmCoordinateX, extWmCoordinateY, extWmCoordinateZ;
    
    // length of the original segment
    length = sqrt( ( pialCoordinateX - wmCoordinateX ) * ( pialCoordinateX - wmCoordinateX ) + ( pialCoordinateY - wmCoordinateY ) * ( pialCoordinateY - wmCoordinateY ) +( pialCoordinateZ - wmCoordinateZ ) * ( pialCoordinateZ - wmCoordinateZ ) );
    
    if ( 0 >= length ) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "*** WARNING: profile length <= zero!!" << std::endl;
        std::cerr << std::endl;
    } else {
        // length of the extention at both sides
        extPartLength = (( length * extentFactor ) - length) / 2;
        
        normX = ( pialCoordinateX - wmCoordinateX ) / length;
        normY = ( pialCoordinateY - wmCoordinateY ) / length;
        normZ = ( pialCoordinateZ - wmCoordinateZ ) / length;
        
        // compute new pial point including extention for pial boundary
        extPialCoordinateX = pialCoordinateX + (normX * extPartLength);
        extPialCoordinateY = pialCoordinateY + (normY * extPartLength);
        extPialCoordinateZ = pialCoordinateZ + (normZ * extPartLength);
        
        // for extention for wm boundary multiply normX etc by -1
        extWmCoordinateX = wmCoordinateX + (-1 * normX * extPartLength);
        extWmCoordinateY = wmCoordinateY + (-1 * normY * extPartLength);
        extWmCoordinateZ = wmCoordinateZ + (-1 * normZ * extPartLength);
        
        // step computed from wm towards pial surface
        // note that we want a certain nr of steps (not segments) so we go from 0 to (nrOfSteps - 1)
        xStep = ( extPialCoordinateX - extWmCoordinateX ) / ( nrOfSteps - 1);
        yStep = ( extPialCoordinateY - extWmCoordinateY ) / ( nrOfSteps - 1);
        zStep = ( extPialCoordinateZ - extWmCoordinateZ ) / ( nrOfSteps - 1);
        
        point[0]=extWmCoordinateX;
        point[1]=extWmCoordinateY;
        point[2]=extWmCoordinateZ;
        
        
        std::cout << lineNr << " " << extPialCoordinateX << " " << extPialCoordinateY << " " << extPialCoordinateZ << " " << extWmCoordinateX << " " << extWmCoordinateY << " " << extWmCoordinateZ << " ";
        
        
        for (i = 0; i < nrOfSteps; i++) {
            
            const bool isInside = inVolume->TransformPhysicalPointToIndex(point, pixelIndex); // test if coordinate falls within image
            
            if ( isInside ) {
                pixelValue = inter->Evaluate(point); // Now use linear interpolation
                std::cout << pixelValue << " ";
            } else {
                std::cerr << std::endl;
                std::cerr << std::endl;
                std::cerr << "*** ERROR: coordinate of sample point is outside volume (should not occur!)" << std::endl;
                std::cerr << std::endl;
                exit(EXIT_FAILURE);
            }
            
            point[0]+=xStep;
            point[1]+=yStep;
            point[2]+=zStep;
        }
        std::cout << length << " " << SCL << std::endl;
    }
}



//
// MAIN LINE
//

int main(int argc, char* argv[]) {
    
    //
    //  Parse program commands with Metacommand
    //
    
    MetaCommand command;
    
    bool useFlagVolume = false;
    bool linearInterpolation = false;
    bool bSplineInterpolation = false;
    bool negateXY = false;
    int bSpline = 2;
    
    command.AddField("wmCoordinates", "name of file containing corresponding white matter boundary coordinates", MetaCommand::STRING,true);
    command.AddField("pialCoordinates", "name of file containing pial coordinates", MetaCommand::STRING,true);
    command.AddField("dataToSample", "name of nifti file with data to sample", MetaCommand::STRING,true);
    
    command.SetOption("nrOfSteps","s",false, "number of steps to sample a profile");
    command.SetOptionLongTag("nrOfSteps","nrOfSteps");
    command.AddOptionField("nrOfSteps", "nos", MetaCommand::INT, true,"100");
    
    command.SetOption("extentFactor","e",false, "increase the segment in both ends with this factor. If 1 (default) then the line will run exactly from layer 1 to layer 2. If larger than 1 the line will extent in both directions (eg 3.0 will give extentions on both ends of the same size as the original line)");
    command.SetOptionLongTag("extentFactor","extentFactor");
    command.AddOptionField("extentFactor", "extent", MetaCommand::INT, true,"0");
    
    command.SetOption("linearInterpolation","l",false, "linear interpolation instead of nearest neighbour");
    command.SetOptionLongTag("linearInterpolation","linearInterpolation");
    
    command.SetOption("bSplineInterpolation","b", false, "b-spline interpolation instead of nearest neighbour (b-spline value 0 - 5)");
    command.SetOptionLongTag("bSplineInterpolation","bSplineInterpolation");
    command.AddOptionField("bSplineInterpolation","bspline", MetaCommand::INT,true,"-100");
    
    command.SetOption("flagVolume","f",false,"write flagvolume for debugging purposes");
    command.SetOptionLongTag("flagVolume","flagVolume");
    command.AddOptionField("flagVolume","flagFile", MetaCommand::STRING,true,"");
    
    command.SetOption("scalarAddOn","a",false, "optional file with scalar info (e.g. curvature) to be added to the measurements");
    command.SetOptionLongTag("scalarAddOn", "scalarAddOn");
    command.AddOptionField("scalarAddOn","addOn", MetaCommand::STRING,true,"");
    
    command.SetOption("negateXY","p",false, "For some reason XY coordinates (e.g. from Freesurfer) need to be multiplied with -1 (appears to be a discrepancy between q/s-form in fileheader and info stored in ITK volume) ");
    command.SetOptionLongTag("negateXY","negateXY");
    
    if (!command.Parse(argc, argv)) {
        std::cerr << std::endl;
        std::cerr << std::endl;
        std::cerr << "With this program you sample cortical profiles based on inner and outer cortical coordinates (two freesufer output files) and a (MRI) scan to be sampled." << std::endl;
        std::cerr << std::endl;
        std::cerr << "format output: | profile nr | XYZ inner | XYZ outer | ... nrOfSteps sample points ... | thickness (Euclidian) | [ scalar addOn ] |" << std::endl;
        std::cerr << std::endl;
        std::cerr << "USAGE:" << std::endl;
        std::cerr << "samplecortex <options> <inner-coordinates>.txt <outer-coordinates>.txt <inputVolume>.nii" << std::endl << std::endl << std::endl;
        return EXIT_FAILURE;
    }
    
    
    std::string pialCoordinatesFileName = command.GetValueAsString("pialCoordinates");
    std::string wmCoordinatesFileName   = command.GetValueAsString("wmCoordinates");
    std::string dataToSampleFileName    = command.GetValueAsString("dataToSample");
    
    int nrOfSteps                       = command.GetValueAsInt("nrOfSteps","nos");
    float extentFactor                  = command.GetValueAsFloat("extentFactor","extent");
    linearInterpolation                 = command.GetValueAsBool("linearInterpolation","linearInterpolation");
    bSpline                             = command.GetValueAsInt("bSplineInterpolation","bspline");
    negateXY                            = command.GetValueAsBool("negateXY","negateXY");
    
    std::string scalarAddOnFileName     = command.GetValueAsString("scalarAddOn", "addOn");
    std::string flagFileName            = command.GetValueAsString("flagVolume","flagFile");
    
    // test if parameter was set
    if (bSpline != -100 ) bSplineInterpolation = true;
    
    // test if both bspline and linear interpolation flags are set
    if ( linearInterpolation && bSplineInterpolation ) {
        std::cerr << std::endl;
        std::cerr << "ERROR: both linear and bSpline interpolation flags are set!" << std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;
    }
    
    if (bSplineInterpolation && (bSpline < 0 || bSpline > 5) ) {
        std::cerr << std::endl;
        std::cerr << "ERROR: bSpline value must be in range 0 - 5" << std::endl;
        std::cerr << std::endl;
        return EXIT_FAILURE;

    }
    
    if (!flagFileName.empty()) useFlagVolume = true;
    
    // EOF command parsing
    
    itk::ImageIOBase::Pointer volumeIO = getImageIO(dataToSampleFileName);
    
    using inVolReader = itk::ImageFileReader<ImageType>;
    FILE* pialCoordinates = NULL;
    FILE* wmCoordinates = NULL;
    FILE* scalarAddOn = NULL;
    
    inVolReader::Pointer in_volReader = inVolReader::New();
    in_volReader->SetFileName(volumeIO->GetFileName());
    in_volReader->Update();
   
    // To create a flagVolume initiated with zeroes we (mis)use a threshold filter.
    using FilterType = itk::BinaryThresholdImageFilter< ImageType, ImageType >;
    FilterType::Pointer filter = FilterType::New();
    
    ImageType::Pointer inVolume;
    inVolume = in_volReader->GetOutput();
    
    
    // Needed if we want to do linear interpolation instead of nearest neighbour
    linearInterpolatorType::Pointer linearInter = linearInterpolatorType::New();
    bSplineInterpolatorType::Pointer bSplineInter = bSplineInterpolatorType::New();
    
    if (linearInterpolation) {
        linearInter->SetInputImage(inVolume);
    }
        
    if (bSplineInterpolation) {
        bSplineInter->SetInputImage(inVolume);
        bSplineInter->SetSplineOrder(bSpline);
    }
    
    ImageType::Pointer flagVolume;
    
    // Note that the flag volume has the coordinate LPS coordinate system if the flag is set
    
    if ( useFlagVolume ) {
        filter->SetInput(inVolume);
        filter->SetOutsideValue(0.0);
        filter->SetInsideValue(0.0);
        filter->SetLowerThreshold(0.0);
        filter->SetUpperThreshold(0.0);
        filter->Update();
        flagVolume = filter->GetOutput();
    }
    //
    // Read in the two ascii files containing the begin and end points of the sample lines, respectively
    //
    
    pialCoordinates=std::fopen(pialCoordinatesFileName.c_str(), "r");
    wmCoordinates=std::fopen(wmCoordinatesFileName.c_str(), "r");
    
    if ( !scalarAddOnFileName.empty() ) {
        scalarAddOn=std::fopen(scalarAddOnFileName.c_str(), "r");
    }
    
    if (std::ferror(pialCoordinates)) {
        std::cout << "file: " << pialCoordinatesFileName << " could not be opened for reading\n\n" << std::endl << std::endl;
        return EXIT_FAILURE;
    }
    if (std::ferror(wmCoordinates)) {
        std::cout << "file: " << wmCoordinatesFileName << " could not be opened for reading\n\n" << std::endl << std::endl;
        return EXIT_FAILURE;
    }

    
    unsigned int i,j,k=1;
    
    float n1,n2,n3;
    float pialCoordinateX, pialCoordinateY, pialCoordinateZ;
    float wmCoordinateX, wmCoordinateY, wmCoordinateZ;
    float SCL;

    long lineNr=1;
    
    while (!feof(pialCoordinates) && !feof(wmCoordinates)) {
        
        i=fscanf(pialCoordinates, "%f %f %f\n", &n1, &n2, &n3);
        pialCoordinateX = n1; pialCoordinateY = n2; pialCoordinateZ = n3;
        
        j=fscanf(wmCoordinates, "%f %f %f\n", &n1, &n2, &n3);
        wmCoordinateX = n1; wmCoordinateY = n2; wmCoordinateZ = n3;
        
        if (negateXY) {
            pialCoordinateX = pialCoordinateX * -1;
            pialCoordinateY = pialCoordinateY * -1;
            wmCoordinateX = wmCoordinateX * -1;
            wmCoordinateY = wmCoordinateY * -1;
        }
        
        
        if (j<3 && i!=j) {
            std::cout << "ERROR!: invalid coordinate (missing value, bailing out)" << std::endl << std::endl;
            return EXIT_FAILURE;
        }

        // do we have scalarAddOns; if so read
        if (scalarAddOn != NULL) {
            k=fscanf(scalarAddOn, "%f\n", &SCL);
        } else SCL = 0;
        
        if (k!=1) {
            std::cout << "ERROR!: invalid scalarAddOn (k = " << k << ") " << std::endl << std::endl;
            return EXIT_FAILURE;
        }

        if (feof(pialCoordinates) != feof(wmCoordinates)) {
            std::cout << "ERROR!!!: pial and wm files differ in length!!" << std::endl << std::endl;
            return EXIT_FAILURE;
        }
        
        if ( linearInterpolation ) {
            if (useFlagVolume) {
                
                sampleAndWriteLineWithFlagsLI(inVolume, linearInter, flagVolume, pialCoordinateX, pialCoordinateY, pialCoordinateZ, wmCoordinateX, wmCoordinateY, wmCoordinateZ, nrOfSteps, extentFactor, lineNr, SCL);
            } else {
                sampleAndWriteLineLI(inVolume, linearInter, pialCoordinateX, pialCoordinateY, pialCoordinateZ, wmCoordinateX, wmCoordinateY, wmCoordinateZ, nrOfSteps, extentFactor, lineNr, SCL);
            }
        } else if ( bSplineInterpolation ) {
            if (useFlagVolume) {
                sampleAndWriteLineWithFlagsBS(inVolume, bSplineInter, flagVolume, pialCoordinateX, pialCoordinateY, pialCoordinateZ, wmCoordinateX, wmCoordinateY, wmCoordinateZ, nrOfSteps, extentFactor, lineNr, SCL);
            } else {
                sampleAndWriteLineBS(inVolume, bSplineInter, pialCoordinateX, pialCoordinateY, pialCoordinateZ, wmCoordinateX, wmCoordinateY, wmCoordinateZ, nrOfSteps, extentFactor, lineNr, SCL);
            }
        } else {
            if (useFlagVolume) {
                sampleAndWriteLineWithFlagsNN(inVolume, flagVolume, pialCoordinateX, pialCoordinateY, pialCoordinateZ, wmCoordinateX, wmCoordinateY, wmCoordinateZ, nrOfSteps, extentFactor, lineNr, SCL );
            } else {
                sampleAndWriteLineNN(inVolume, pialCoordinateX, pialCoordinateY, pialCoordinateZ, wmCoordinateX, wmCoordinateY, wmCoordinateZ, nrOfSteps, extentFactor, lineNr, SCL );
            }
        }
        lineNr++;
    }
    
     std::string outName = flagFileName;
    
    // Write out the flag volume if any
    if (useFlagVolume) {
        using WriterType = itk::ImageFileWriter<ImageType>;
        try {
            WriterType::Pointer writer = WriterType::New();
            writer->SetFileName(outName);
            writer->SetInput(flagVolume);
            writer->Update();
        } catch ( itk::ExceptionObject & e ) {
            std::cerr << "Unexpected exception caught when writing. Note, requires valid file extention (e.g. .nii): " << outName << std::endl;
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
