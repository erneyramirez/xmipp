/***************************************************************************
 *
 * Authors:     Jeison Méndez García (jmendez@utp.edu.co)
 *
 * Instituto de Investigaciones en Matemáticas Aplicadas y en Sistemas -- IIMAS
 * Universidad Nacional Autónoma de México -UNAM
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "angular_assignment_mag.h"

ProgAngularAssignmentMag::ProgAngularAssignmentMag()
{
    produces_a_metadata = true;
    each_image_produces_an_output = false;
}

ProgAngularAssignmentMag::~ProgAngularAssignmentMag()
{

}

void ProgAngularAssignmentMag::defineParams()
{
    XmippMetadataProgram::defineParams();
    //usage
    addUsageLine("Generates a list of candidates for angular assignment for each experimental image");
    //params
//    addParamsLine("   -i <md_file>               : Metadata file with input experimental projections");
//    addParamsLine("   -o <md_file>               : Metadata file with output projections");
    addParamsLine("   -ref <md_file>             : Metadata file with input reference projections");
    addParamsLine("  [-odir <outputDir=\".\">]   : Output directory");
    addParamsLine("  [-sym <symfile=c1>]         : Enforce symmetry in projections");
    addParamsLine("  [-sampling <sampling=1.>]         : sampling");
}

// Read arguments ==========================================================
void ProgAngularAssignmentMag::readParams()
{
    XmippMetadataProgram::readParams();
    fnIn = XmippMetadataProgram::fn_in;
    fnOut = XmippMetadataProgram::fn_out;
    fnRef = getParam("-ref");
    fnDir = getParam("-odir");
    sampling = getDoubleParam("-sampling");
    XmippMetadataProgram::oroot = fnDir;
    fnSym = getParam("-sym");
}

// Show ====================================================================
void ProgAngularAssignmentMag::show()
{
    if (verbose > 0)
    {
	printf("%d reference images of %d x %d\n", int(sizeMdRef), int(Xdim), int(Ydim));
        printf("%d exp images of %d x %d in this group\n", int(sizeMdIn), int(Xdim), int(Ydim));
        //        printf("imgcc %d x %d from mdIn:%d, mdRef:%d\n", int(YSIZE(imgcc)), int(XSIZE(imgcc)), int(sizeMdIn), int(sizeMdRef));
        printf("\nstartBand= %d\n", int(startBand));
        printf("finalBand= %d\n", int(finalBand3));
        printf("n_bands= %d\n", int(n_bands3));

	XmippMetadataProgram::show();
	//        std::cout << "Input metadata              : "  << fnIn        << std::endl;
	std::cout << "Input references: "  << fnRef       << std::endl;
	//        std::cout << "Output directory            : "  << fnDir       << std::endl;
	//        if (fnSym != "")
	//            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
	std::cout << "sampling: " << sampling << std::endl;
    }
}

void ProgAngularAssignmentMag::startProcessing()
{
    XmippMetadataProgram::startProcessing();
}

/* print in console some values of double MultidimArray */
void ProgAngularAssignmentMag::printSomeValues(MultidimArray<double> &MDa){
    for(int i=0; i<4; i++)
	for(int j=0; j<4; j++)
	    std::cout << "val: " << DIRECT_A2D_ELEM(MDa,i,j) << std::endl;
}

void ProgAngularAssignmentMag::preProcess()
{
    mdIn.read(fnIn);
    mdRef.read(fnRef);

    //borrar esta
    cont=0;

    // size of images
    size_t Zdim, Ndim;
    getImageSize(mdIn,Xdim,Ydim,Zdim,Ndim);

    // some constants
    n_rad = size_t(Xdim/2 + 0.5);

    // estos indices de inicio y parada van a ser los de los extremos
    //tercer conjunto de bandas
    startBand=size_t((sampling*Xdim)/40.);
    finalBand3=size_t((sampling*Xdim)/16.); // este valor debe estar relacionado con el maxResolutionTarget del protocolo
    finalBand3=(finalBand3 >= n_rad) ? n_rad-1 : finalBand3;
    n_bands3=finalBand3-startBand;

    // primer conjunto de bandas
    n_bands1=int(n_bands3/3.+.5);
    finalBand1=startBand+n_bands1;

    // segundo conjunto de bandas
    n_bands2=2*n_bands1;
    finalBand2=startBand+n_bands2;

    // non-overlapped bands
    startBand2=finalBand2-n_bands1;
    startBand3=finalBand3-n_bands1;
    n_bands3=finalBand3-startBand3;
    n_bands2=finalBand2-startBand2;

    //    // overlapped bands
    //    startBand2=startBand;
    //    startBand3=startBand;

    //definir para operar bandas por separado -- No traslapar


    printf("start-final:[%d-%d] y ancho %d\n",int(startBand),int(finalBand1), int(n_bands1));
    printf("start-final:[%d-%d] y ancho %d\n",int(startBand2),int(finalBand2), int(n_bands2));
    printf("start-final:[%d-%d] y ancho %d\n",int(startBand3),int(finalBand3), int(n_bands3));


    //angulos
    n_ang = size_t(180);
    n_ang2 = 2*n_ang;
    maxShift = .10 * Xdim; // read maxShift as input parameter

    // read reference images
    FileName fnImgRef;
    MDRow rowRef;
    sizeMdRef = mdRef.size();

    // how many input images
    sizeMdIn = mdIn.size();

    // reference image related
    Image<double>                           ImgRef;
    MultidimArray<double>                   MDaRef(Ydim,Xdim);
    MultidimArray< std::complex<double> >   MDaRefF  ;
    MultidimArray< std::complex<double> >   MDaRefF2 ;
    MultidimArray<double>                   MDaRefFM ;
    MultidimArray<double>                   MDaRefFMs;
    MultidimArray<double>                   MDaRefFMs_polarPart1(n_bands1, n_ang2);
    MultidimArray< std::complex<double> >   MDaRefFMs_polarF1;
    MultidimArray<double>                   MDaRefFMs_polarPart2(n_bands2, n_ang2);
    MultidimArray< std::complex<double> >   MDaRefFMs_polarF2;
    MultidimArray<double>                   MDaRefFMs_polarPart3(n_bands3, n_ang2);
    MultidimArray< std::complex<double> >   MDaRefFMs_polarF3;

    // try to storage all data related to reference images in memory
    for (int k = 0; k < sizeMdRef; k++){
	// reading image
	mdRef.getRow(rowRef, size_t(k+1) );
	rowRef.getValue(MDL_IMAGE, fnImgRef);
	// processing reference image
	ImgRef.read(fnImgRef);
	MDaRef = ImgRef();
	vecMDaRef.push_back(MDaRef);
	_applyFourierImage2(MDaRef, MDaRefF);
	vecMDaRefF.push_back(MDaRefF);
	transformerImage.getCompleteFourier(MDaRefF2);
	_getComplexMagnitude(MDaRefF2, MDaRefFM);
        completeFourierShift(MDaRefFM, MDaRefFMs);
        //aplica fourier y almacena en el ancho de banda 1
        MDaRefFMs_polarPart1 = imToPolar(MDaRefFMs,startBand,finalBand1,n_bands1);
        _applyFourierImage2(MDaRefFMs_polarPart1, MDaRefFMs_polarF1, n_ang);
        vecMDaRefFMs_polarF1.push_back(MDaRefFMs_polarF1);
        //        printf("size MDaRefFMs_polarF1: %d-%d\n",int(YSIZE(MDaRefFMs_polarF1)),int(XSIZE(MDaRefFMs_polarF1)));
        //aplica fourier y almacena en el ancho de banda 2
        MDaRefFMs_polarPart2 = imToPolar(MDaRefFMs,startBand2,finalBand2,n_bands2);
        _applyFourierImage2(MDaRefFMs_polarPart2, MDaRefFMs_polarF2, n_ang);
        vecMDaRefFMs_polarF2.push_back(MDaRefFMs_polarF2);
        //        printf("size MDaRefFMs_polarF2: %d-%d\n",int(YSIZE(MDaRefFMs_polarF2)),int(XSIZE(MDaRefFMs_polarF2)));
        //aplica fourier y almacena en el ancho de banda 3
        MDaRefFMs_polarPart3 = imToPolar(MDaRefFMs,startBand3,finalBand3,n_bands3);
        _applyFourierImage2(MDaRefFMs_polarPart3, MDaRefFMs_polarF3, n_ang);
        vecMDaRefFMs_polarF3.push_back(MDaRefFMs_polarF3);
        //        printf("size MDaRefFMs_polarF3: %d-%d\n",int(YSIZE(MDaRefFMs_polarF3)),int(XSIZE(MDaRefFMs_polarF3)));
    }

    //    printf("terminó con imágenes referencia\n");
    //    exit(1);

    // // creo que debo definir esto luego dentro de processImage() quizá también sea problematico??
    // aunque cuando pongo printf y MPI entonces imprime varias veces... quizá no esté mal ponerlos acá
    candidatesFirstLoop.resize(sizeMdRef);
    Idx.resize(sizeMdRef);
    candidatesFirstLoopCoeff.resize(sizeMdRef);
    bestTx.resize(sizeMdRef);
    bestTy.resize(sizeMdRef);
    bestRot.resize(sizeMdRef);

    candidatesSecondLoop.resize(sizeMdRef);
    Idx2.resize(sizeMdRef);
    candidatesSecondLoopCoeff.resize(sizeMdRef);
    bestTx2.resize(sizeMdRef);
    bestTy2.resize(sizeMdRef);
    bestRot2.resize(sizeMdRef);

    candidatesThirdLoop.resize(sizeMdRef);
    Idx3.resize(sizeMdRef);
    candidatesThirdLoopCoeff.resize(sizeMdRef);
    bestTx3.resize(sizeMdRef);
    bestTy3.resize(sizeMdRef);
    bestRot3.resize(sizeMdRef);
    // delay axes
    //    _delayAxes(Ydim, Xdim, n_ang);
    mdOut.setComment("experiment for metadata output containing data for reconstruction");
}

void ProgAngularAssignmentMag::processImage(const FileName &fnImg, const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut){

    // experimental image related
    rowOut = rowIn;

    // input image related
    MDRow rowRef;
    Image<double>                           ImgIn;
    MultidimArray<double>                   MDaIn(Ydim,Xdim);
    MultidimArray< std::complex<double> >   MDaInF  ;
    MultidimArray< std::complex<double> >   MDaInF2 ;
    MultidimArray<double>                   MDaInFM ;
    MultidimArray<double>                   MDaInFMs;
    //n_bands1
    MultidimArray<double>                   MDaInFMs_polarPart1(n_bands1, n_ang2);
    MultidimArray< std::complex<double> >   MDaInFMs_polarF1; // mirar luego si esta se puede reutilizar
    //n_bands2
    MultidimArray<double>                   MDaInFMs_polarPart2(n_bands2, n_ang2);
    MultidimArray< std::complex<double> >   MDaInFMs_polarF2; // mirar luego si esta se puede reutilizar
    //n_bands3
    MultidimArray<double>                   MDaInFMs_polarPart3(n_bands3, n_ang2);
    MultidimArray< std::complex<double> >   MDaInFMs_polarF3; // mirar luego si esta se puede reutilizar

    // processing input image
    ImgIn.read(fnImg);
    MDaIn = ImgIn();
    _applyFourierImage2(MDaIn, MDaInF);
    transformerImage.getCompleteFourier(MDaInF2);
    _getComplexMagnitude(MDaInF2, MDaInFM);
    completeFourierShift(MDaInFM, MDaInFMs);

    // BANDA 1
    // a partir de aquí es diferente para cada resolución.
    MDaInFMs_polarPart1 = imToPolar(MDaInFMs,startBand,finalBand1,n_bands1);
    _applyFourierImage2(MDaInFMs_polarPart1, MDaInFMs_polarF1, n_ang);

    tempCoeff = 0.0;
    int k = 0;
    double bestCandVar, bestCoeff, Tx, Ty;

    std::ofstream outfile("/home/jeison/Escritorio/TestOutfile.txt");

    double band1_bestCoeff=-1.0;
    int idx_band1=0;

    // loop over reference stack
    for(int countRefImg = 0; countRefImg < sizeMdRef; countRefImg++){
        // computing relative rotation and traslation
        ccMatrix(MDaInFMs_polarF1, vecMDaRefFMs_polarF1[countRefImg], ccMatrixRot);
        maxByColumn(ccMatrixRot, ccVectorRot); // el meanByColumn ha dado segmentation fault en asimov
        peaksFound = 0;
        std::vector<double>().swap(cand);
        rotCandidates3(ccVectorRot, cand, XSIZE(ccMatrixRot), &peaksFound); // rotcandidates3 // aquí hay que revisar lo que puedo hacer cuando peaksFound==0
        bestCand(MDaIn, MDaInF, vecMDaRef[countRefImg], cand, peaksFound, &bestCandVar, &Tx, &Ty, &bestCoeff);
        // all the results are storaged for posterior partial_sort
        Idx[countRefImg] = k++;
        candidatesFirstLoop[countRefImg] = countRefImg+1;
        candidatesFirstLoopCoeff[countRefImg] = bestCoeff;
        bestTx[countRefImg] = Tx;
        bestTy[countRefImg] = Ty;
        bestRot[countRefImg] = bestCandVar;
        if(bestCoeff>band1_bestCoeff){
            band1_bestCoeff=bestCoeff;
            idx_band1=Idx[countRefImg];
        }
    }



    // BANDA 2
    // a partir de aquí es diferente para cada resolución
    MDaInFMs_polarPart2 = imToPolar(MDaInFMs,startBand2,finalBand2,n_bands2);
    _applyFourierImage2(MDaInFMs_polarPart2, MDaInFMs_polarF2, n_ang);

    tempCoeff = 0.0;
    k = 0;
    //bestCandVar, bestCoeff, Tx, Ty;
    // loop over reference stack
    double band2_bestCoeff=-1.0;
    int idx_band2=0;
    for(int countRefImg = 0; countRefImg < sizeMdRef; countRefImg++){
        // computing relative rotation and traslation
        ccMatrix(MDaInFMs_polarF2, vecMDaRefFMs_polarF2[countRefImg], ccMatrixRot);
        maxByColumn(ccMatrixRot, ccVectorRot); // el meanByColumn ha dado segmentation fault (según yo, la suma es muy alta, hay que multiplicar por 1/size cada vez)
        peaksFound = 0;
        std::vector<double>().swap(cand);
        rotCandidates3(ccVectorRot, cand, XSIZE(ccMatrixRot), &peaksFound); // rotcandidates3 // aquí hay que revisar lo que puedo hacer cuando peaksFound==0
        bestCand(MDaIn, MDaInF, vecMDaRef[countRefImg], cand, peaksFound, &bestCandVar, &Tx, &Ty, &bestCoeff);
        // all the results are storaged for posterior partial_sort
        Idx2[countRefImg] = k++;
        candidatesSecondLoop[countRefImg] = countRefImg+1;
        candidatesSecondLoopCoeff[countRefImg] = bestCoeff;
        bestTx2[countRefImg] = Tx;
        bestTy2[countRefImg] = Ty;
        bestRot2[countRefImg] = bestCandVar;
        if(bestCoeff>band2_bestCoeff){
            band2_bestCoeff=bestCoeff;
            idx_band2=Idx2[countRefImg];
        }
    }



    // BANDA 3
    // a partir de aquí es diferente para cada resolución
    MDaInFMs_polarPart3 = imToPolar(MDaInFMs,startBand3,finalBand3,n_bands3);
    _applyFourierImage2(MDaInFMs_polarPart3, MDaInFMs_polarF3, n_ang);

    tempCoeff = 0.0;
    k = 0;
    //bestCandVar, bestCoeff, Tx, Ty;
    // loop over reference stack
    double band3_bestCoeff=-1.0;
    int idx_band3=0;
    for(int countRefImg = 0; countRefImg < sizeMdRef; countRefImg++){
        // computing relative rotation and traslation
        ccMatrix(MDaInFMs_polarF3, vecMDaRefFMs_polarF3[countRefImg], ccMatrixRot);
        maxByColumn(ccMatrixRot, ccVectorRot); // el meanByColumn ha dado segmentation fault en asimov
        peaksFound = 0;
        std::vector<double>().swap(cand);
        rotCandidates3(ccVectorRot, cand, XSIZE(ccMatrixRot), &peaksFound); // rotcandidates3 // aquí hay que revisar lo que puedo hacer cuando peaksFound==0
        bestCand(MDaIn, MDaInF, vecMDaRef[countRefImg], cand, peaksFound, &bestCandVar, &Tx, &Ty, &bestCoeff);
        // all the results are storaged for posterior partial_sort
        Idx3[countRefImg] = k++;
        candidatesThirdLoop[countRefImg] = countRefImg+1;
        candidatesThirdLoopCoeff[countRefImg] = bestCoeff;
        bestTx3[countRefImg] = Tx;
        bestTy3[countRefImg] = Ty;
        bestRot3[countRefImg] = bestCandVar;
        if(bestCoeff>band3_bestCoeff){
            band3_bestCoeff=bestCoeff;
            idx_band3=Idx3[countRefImg];
        }
    }

    // producto de la correlación normalizada para las 3 bandas procesadas
    // Idx1-2-3 son los mismos. se puede aquí crear uno solo de esos
    std::vector<double> candidatesNormalizedProduct;
    candidatesNormalizedProduct.resize(sizeMdRef);
    for(int ii=0;ii<sizeMdRef;ii++){
        candidatesNormalizedProduct[ii]=(candidatesFirstLoopCoeff[ii]/band1_bestCoeff) *
                                        (candidatesSecondLoopCoeff[ii]/band2_bestCoeff) *
                                        (candidatesThirdLoopCoeff[ii]/band3_bestCoeff);
    }

    int nCand = sizeMdRef; // 1  3
    std::partial_sort(Idx.begin(), Idx.begin()+nCand, Idx.end(),
                      [&](int i, int j){return candidatesNormalizedProduct[i] > candidatesNormalizedProduct[j]; });

    /*
    int howMany=6;
    if(cont<=howMany){
        printf("maxBand1-idx: %.4f-%d,  maxBand2-idx: %.4f-%d,   maxBand3-idx: %.4f-%d, \n",
               band1_bestCoeff, idx_band1,
               band2_bestCoeff, idx_band2,
               band3_bestCoeff, idx_band3);
        printf("en candidato %d \t band1_corr %.4f  \t band2_corr %.4f  \t band3_corr %.4f\n",
                candidatesFirstLoop[Idx[0]],
                candidatesFirstLoopCoeff[Idx[0]],
                candidatesSecondLoopCoeff[Idx[0]],
                candidatesThirdLoopCoeff[Idx[0]]);
        printf("en candidato %d \t band1_Rot %.2f  \t band2_Rot %.2f  \t band3_Rot %.2f\n",
                candidatesFirstLoop[Idx[0]],
                bestRot[Idx[0]],
                bestRot2[Idx[0]],
                bestRot3[Idx[0]]);
        printf("en candidato %d \t band1_tx_ty [%.2f, %.2f]  \t band2_tx_ty [%.2f, %.2f]  \t band3_tx_ty [%.2f,%.2f]\n",
                candidatesFirstLoop[Idx[0]],
                bestTx[Idx[0]], bestTy[Idx[0]],
                bestTx2[Idx[0]], bestTy2[Idx[0]],
                bestTx3[Idx[0]], bestTy3[Idx[0]]);
    }
    if (cont==howMany){
        // test de escribir los resultados
        for(int i=0;i<nCand;i++){
            outfile<<candidatesFirstLoopCoeff[i]
                   <<"\t"<<candidatesSecondLoopCoeff[i]
                   <<"\t"<<candidatesThirdLoopCoeff[i]
                   <<"\t"<<candidatesNormalizedProduct[i]
                   <<"\t"<<bestRot[i]<<"\t"<<bestRot2[i]<<"\t"<<bestRot3[i]<<"\n\n";
        }
        outfile.close();

        //          // no ha sido muy util porque ahora mismo trabajo solo con experimentales
        //        // imprimir en archivos las imágenes implicadas.
        //        printf("almacena Referencia-mejor_candidato. Idx:%d, rot: %.2f, Tx: %.2f, Tx: %.2f\n",
        //               candidatesFirstLoop[Idx[0]]-1,bestRot[Idx[0]],bestTx[Idx[0]], bestTy[Idx[0]]);
        //        // aplicar transformación e imprimir
        //        MultidimArray<double> MDaTransformada;
        //        _applyRotationAndShift(vecMDaRef[candidatesFirstLoop[Idx[0]-1]],bestRot[Idx[0]],
        //                               bestTx[Idx[0]], bestTy[Idx[0]],MDaTransformada);
        //        _writeTestFile(MDaTransformada,"/home/jeison/Escritorio/t_referenciaTransformada.txt",
        //                      YSIZE(vecMDaRef[candidatesFirstLoop[Idx[0]-1]]),
        //                      XSIZE(vecMDaRef[candidatesFirstLoop[Idx[0]-1]]));

        //        printf("almacena experimental\n");
        //        _writeTestFile(MDaIn,"/home/jeison/Escritorio/t_experimental.txt",
        //                      YSIZE(MDaIn),
        //                      XSIZE(MDaIn));

        // solo para la imagen referencia encontrada como la mejor,
        // volver a buscar rotación y shift
        int bestIdxRef=candidatesFirstLoop[Idx[0]]-1;
        MultidimArray<double>  finalRef(Ydim,Xdim); //UNUSED
        MultidimArray< std::complex<double> >   finalRefF  ;
        MultidimArray< std::complex<double> >   finalRefF2 ;
        MultidimArray<double>                   finalRefFM ;
        MultidimArray<double>                   finalRefFMs;
        //n_bands
        size_t finalNbands= finalBand3-startBand;
        printf("finalNbands: %d\n",int(finalNbands));
        MultidimArray<double>                   finalRefFMs_polarPart(finalNbands, n_ang2);
        MultidimArray< std::complex<double> >   finalRefFMs_polarF;

        _applyFourierImage2(vecMDaRef[candidatesFirstLoop[bestIdxRef]], finalRefF);
        transformerImage.getCompleteFourier(finalRefF2);
        _getComplexMagnitude(finalRefF2, finalRefFM);
        completeFourierShift(finalRefFM, finalRefFMs);


        finalRefFMs_polarPart = imToPolar(finalRefFMs,startBand,finalBand3,finalNbands);
        _writeTestFile(finalRefFMs_polarPart,"/home/jeison/Escritorio/z_polarRef.txt",YSIZE(finalRefFMs_polarPart),XSIZE(finalRefFMs_polarPart));
        _applyFourierImage2(finalRefFMs_polarPart, finalRefFMs_polarF, n_ang);

        //experimental
        MultidimArray<double>   MDaInFMs_finalpolarPart(finalNbands, n_ang2);
        MultidimArray< std::complex<double> >   MDaInfinalFMs_polarF;
        MDaInFMs_finalpolarPart = imToPolar(MDaInFMs,startBand,finalBand3,finalNbands);
        _writeTestFile(MDaInFMs_finalpolarPart,"/home/jeison/Escritorio/z_polarEntrada.txt",YSIZE(MDaInFMs_finalpolarPart),XSIZE(MDaInFMs_finalpolarPart));
        _applyFourierImage2(MDaInFMs_finalpolarPart, MDaInfinalFMs_polarF, n_ang);

        ccMatrix(MDaInfinalFMs_polarF, finalRefFMs_polarF, ccMatrixRot);
        maxByColumn(ccMatrixRot, ccVectorRot); // el meanByColumn ha dado segmentation fault en asimov

        std::vector<double>().swap(cand);
        peaksFound=0; bestCandVar=0;
        rotCandidates3(ccVectorRot, cand, XSIZE(ccMatrixRot), &peaksFound); // rotcandidates3 // aquí hay que revisar lo que puedo hacer cuando peaksFound==0
        _writeTestFile(ccVectorRot,"/home/jeison/Escritorio/z_ccVector.txt",1,XSIZE(ccVectorRot));
        bestCand(MDaIn, MDaInF, vecMDaRef[candidatesFirstLoop[bestIdxRef]], cand, peaksFound, &bestCandVar, &Tx, &Ty, &bestCoeff);

        //results
        printf("bestIdxRef+1: %d bestCoeff: %.3f Tx: %.3f Ty: %.3f bestCandVar: %.3f\n",bestIdxRef+1, bestCoeff, Tx, Ty, bestCandVar);


        printf("\n\n*******************detener**************************\n\n");
        exit(1);
    } // */

    // solo para la imagen referencia encontrada como la mejor,
    // volver a buscar rotación y shift
    int bestIdxRef=candidatesFirstLoop[Idx[0]]-1; // Idx fue organizada de mayor a menor basada en valor de correlación
    MultidimArray< std::complex<double> >   finalRefF  ;
    MultidimArray< std::complex<double> >   finalRefF2 ;
    MultidimArray<double>                   finalRefFM ;
    MultidimArray<double>                   finalRefFMs;
    //n_bands
    size_t finalNbands= finalBand3-startBand;
    MultidimArray<double>                   finalRefFMs_polarPart(finalNbands, n_ang2);
    MultidimArray< std::complex<double> >   finalRefFMs_polarF;

    _applyFourierImage2(vecMDaRef[candidatesFirstLoop[bestIdxRef]], finalRefF);
    transformerImage.getCompleteFourier(finalRefF2);
    _getComplexMagnitude(finalRefF2, finalRefFM);
    completeFourierShift(finalRefFM, finalRefFMs);
    finalRefFMs_polarPart = imToPolar(finalRefFMs,startBand,finalBand3,finalNbands);
    _applyFourierImage2(finalRefFMs_polarPart, finalRefFMs_polarF, n_ang);

    //experimental
    MultidimArray<double>   MDaInFMs_finalpolarPart(finalNbands, n_ang2);
    MultidimArray< std::complex<double> >   MDaInfinalFMs_polarF;
    MDaInFMs_finalpolarPart = imToPolar(MDaInFMs,startBand,finalBand3,finalNbands);
    _applyFourierImage2(MDaInFMs_finalpolarPart, MDaInfinalFMs_polarF, n_ang);

    ccMatrix(MDaInfinalFMs_polarF, finalRefFMs_polarF, ccMatrixRot);
    maxByColumn(ccMatrixRot, ccVectorRot); // el meanByColumn ha dado segmentation fault en asimov

    std::vector<double>().swap(cand);
    peaksFound=0; bestCandVar=0;
    rotCandidates3(ccVectorRot, cand, XSIZE(ccMatrixRot), &peaksFound); // rotcandidates3 // aquí hay que revisar lo que puedo hacer cuando peaksFound==0
    bestCand(MDaIn, MDaInF, vecMDaRef[candidatesFirstLoop[bestIdxRef]], cand, peaksFound, &bestCandVar, &Tx, &Ty, &bestCoeff);

    //results
    printf("bestIdxRef+1: %d bestCoeff: %.3f Tx: %.3f Ty: %.3f bestCandVar: %.3f\n",bestIdxRef+1, bestCoeff, Tx, Ty, bestCandVar);



    //        /*
    // skip second loop
   //     choose nCand of the candidates with best corrCoeff
    /*int*/nCand = 1; // 1  3

    double rotRef, tiltRef;
    for(int i = 0; i < nCand; i++){
        // reading info of reference image candidate
        mdRef.getRow(rowRef, size_t( bestIdxRef+1 ) );
        rowRef.getValue(MDL_ANGLE_ROT, rotRef);
        rowRef.getValue(MDL_ANGLE_TILT, tiltRef);

        //save metadata of images with angles
        rowOut.setValue(MDL_IMAGE,       fnImgOut);
        rowOut.setValue(MDL_ENABLED,     1);
        rowOut.setValue(MDL_IDX,         size_t(bestIdxRef+1));
        rowOut.setValue(MDL_MAXCC,       bestCoeff);
        rowOut.setValue(MDL_WEIGHT,      1.);
        rowOut.setValue(MDL_WEIGHT_SIGNIFICANT,   1.);
        rowOut.setValue(MDL_ANGLE_ROT,   rotRef);
        rowOut.setValue(MDL_ANGLE_TILT,  tiltRef);
        rowOut.setValue(MDL_ANGLE_PSI,   bestCandVar);
        rowOut.setValue(MDL_SHIFT_X,     -1. * Tx);
        rowOut.setValue(MDL_SHIFT_Y,     -1. * Ty);
    }
    // */


    //    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    //    std::cout << "Operation took "<< duration*1000 << "milliseconds" << std::endl;

    //borrar
    cont++;
    printf("cont: %d\n",cont);
}

void ProgAngularAssignmentMag::postProcess(){

    // from angularContinousAssign2
    MetaData &ptrMdOut=*getOutputMd();
    ptrMdOut.removeDisabled();
    double maxCC=-1.;
    FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
    {
        double thisMaxCC;
        ptrMdOut.getValue(MDL_MAXCC,thisMaxCC,__iter.objId);
        if (thisMaxCC>maxCC)
            maxCC=thisMaxCC;
    }
    FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
    {
        double thisMaxCC;
        ptrMdOut.getValue(MDL_MAXCC,thisMaxCC,__iter.objId);
        ptrMdOut.setValue(MDL_WEIGHT,thisMaxCC/maxCC,__iter.objId);
        ptrMdOut.setValue(MDL_WEIGHT_SIGNIFICANT,thisMaxCC/maxCC,__iter.objId);
    }

    ptrMdOut.write(XmippMetadataProgram::fn_out.replaceExtension("xmd"));
    transformerImage.cleanup();
    transformerPolarImage.cleanup();
}

/* Pearson Coeff */
void ProgAngularAssignmentMag::pearsonCorr(MultidimArray<double> &X, MultidimArray<double> &Y, double &coeff){

    //    MultidimArray<double>   X2(Ydim,Xdim);
    //    MultidimArray<double>   Y2(Ydim,Xdim);
    //    _applyCircularMask(X,X2);
    //    _applyCircularMask(Y,Y2);
    // covariance
    double X_m, Y_m, X_std, Y_std;
    arithmetic_mean_and_stddev(X, X_m, X_std);
    arithmetic_mean_and_stddev(Y, Y_m, Y_std);

    double prod_mean = mean_of_products(X, Y);
    double covariace = prod_mean - (X_m * Y_m);

    coeff = covariace / (X_std * Y_std);
}

void ProgAngularAssignmentMag::_applyCircularMask(const MultidimArray<double> &in, MultidimArray<double> &out){

    double Cf = (Ydim + (Ydim % 2)) / 2.0; // for odd/even images
    double Cc = (Xdim + (Xdim % 2)) / 2.0;
    int pixReduc = 1;
    double rad2 = (Cf - pixReduc) * (Cf - pixReduc);
    double val = 0;
    out.initZeros(Ydim,Xdim);
    for(size_t f = 0; f < Ydim; f++){
	for(size_t c = 0; c < Xdim; c++){
	    val = (f-Cf)*(f-Cf) + (c-Cc)*(c-Cc);
	    if (val < rad2)
		DIRECT_A2D_ELEM(out, f, c) = DIRECT_A2D_ELEM(in,f,c);
	}
    }
}

/* Arithmetic mean and stdDev for Pearson Coeff */
void ProgAngularAssignmentMag::arithmetic_mean_and_stddev( MultidimArray<double> &data, double &avg, double &stddev ){
    data.computeAvgStdev(avg, stddev);
}

/* Mean of products for Pearson Coeff */
double ProgAngularAssignmentMag::mean_of_products(MultidimArray<double> &data1, MultidimArray<double> &data2){
    double total = 0;
    for (int f = 0; f < Ydim; f++){
	for (int c = 0; c < Xdim; c++){
	    total += DIRECT_A2D_ELEM(data1,f,c) * DIRECT_A2D_ELEM(data2,f,c);
	}
    }
    return total/(Xdim*Ydim);
}

/* writing out some data to file with an specified size*/
void ProgAngularAssignmentMag::_writeTestFile(MultidimArray<double> &data, const char* fileName,
					      size_t nFil, size_t nCol){
    std::ofstream outFile(fileName);
    for (int f = 0; f < nFil; f++){
	for (int c = 0; c < nCol; c++){
	    outFile <<  DIRECT_A2D_ELEM(data,f,c) << "\t";
	}
	outFile << "\n";
    }
    outFile.close();
}

/* writing out some data to file Ydim x Xdim size*/
void ProgAngularAssignmentMag::_writeTestFile(MultidimArray<double> &data, const char* fileName){
    std::ofstream outFile(fileName);
    for (int f = 0; f < Ydim; f++){
	for (int c = 0; c < Xdim; c++){
	    outFile <<  DIRECT_A2D_ELEM(data,f,c) << "\t";
	}
	outFile << "\n";
    }
    outFile.close();
}

/* get COMPLETE fourier spectrum of Images. It should be changed for half */
void ProgAngularAssignmentMag::_applyFourierImage(MultidimArray<double> &data,
						  MultidimArray< std::complex<double> > &FourierData){
    transformerImage.completeFourierTransform(data, FourierData);
}

/* get COMPLETE fourier spectrum of polarRepresentation of Magnitude. It should be changed for half */
void ProgAngularAssignmentMag::_applyFourierImage(MultidimArray<double> &data,
						  MultidimArray< std::complex<double> > &FourierData, const size_t &ang){
    transformerPolarImage.completeFourierTransform(data, FourierData);
}

/*first try in using only one half of Fourier space*/
void ProgAngularAssignmentMag::_applyFourierImage2(MultidimArray<double> &data,
						   MultidimArray< std::complex<double> > &FourierData){
    transformerImage.FourierTransform(data,FourierData,true);
}

/* first try one half of fourier spectrum of polarRepresentation of Magnitude*/
void ProgAngularAssignmentMag::_applyFourierImage2(MultidimArray<double> &data,
						  MultidimArray< std::complex<double> > &FourierData, const size_t &ang){
    transformerPolarImage.FourierTransform(data,FourierData,true); // false --> true para generar copia
}


/* get magnitude of fourier spectrum */
void ProgAngularAssignmentMag::_getComplexMagnitude( MultidimArray< std::complex<double> > &FourierData,
						     MultidimArray<double> &FourierMag){
    FFT_magnitude(FourierData,FourierMag);
}

/* cartImg contains cartessian  grid representation of image,
*  rad and ang are the number of radius and angular elements*/
MultidimArray<double> ProgAngularAssignmentMag::imToPolar(MultidimArray<double> &cartIm,
							  size_t &startBand,
                                                          size_t &finalBand,
                                                          size_t &n_bands){

    MultidimArray<double> polarImg(n_bands, n_ang2);
    float pi = 3.141592653;
    // coordinates of center
    //    double cy = (Ydim+1)/2.0;
    //    double cx = (Xdim+1)/2.0;
    double cy = (Ydim + (Ydim % 2)) / 2.0; // for odd/even images
    double cx = (Xdim + (Xdim % 2)) / 2.0;

    // scale factors
    double sfy = (Ydim-1)/2.0;
    double sfx = (Xdim-1)/2.0;

    double delR = (double)(1.0 / (n_rad)); // n_rad-1
    double delT = 2.0 * pi / n_ang2;

    // loop through rad and ang coordinates
    double r, t, x_coord, y_coord;
    for(size_t ri = startBand; ri < finalBand; ri++){
	for(size_t ti = 0; ti < n_ang2; ti++ ){
	    r = ri * delR;
	    t = ti * delT;
	    x_coord = ( r * cos(t) ) * sfx + cx;
	    y_coord = ( r * sin(t) ) * sfy + cy;
	    // set value of polar img
	    DIRECT_A2D_ELEM(polarImg,ri-startBand,ti) = interpolate(cartIm,x_coord,y_coord);
	}
    }

//    printf("termina polar\n r, t, xcoord, ycoord = %.2f, %.2f, %.2f, %.2f\n", r, t, x_coord, y_coord);

    return polarImg;
}

/* cartImg contains cartessian  grid representation of image,
*  rad and ang are the number of radius and angular elements
*  this function was built for half representation of Fourier spectrum*/
MultidimArray<double> ProgAngularAssignmentMag::imToPolar2(MultidimArray<double> &cartIm,
							  const size_t &rad, const size_t &ang){
    MultidimArray<double> polarImg(rad, ang);
    float pi = 3.141592653;
    // coordinates of center
    double cy = 0.5; //(Ydim+1)/2.0;
    double cx = (Xdim+1)/2.0;
    // scale factors
    double sfy = (Ydim-1)/2.0;
    double sfx = (Xdim-1)/2.0;

    double delR = (double)(1.0 / (rad-1));
    double delT = pi / ang;

    // loop through rad and ang coordinates
    double r, t, x_coord, y_coord;
    for(size_t ri = 0; ri < rad; ri++){
	for(size_t ti = 0; ti < ang; ti++ ){
	    r = ri * delR;
	    t = ti * delT;
	    x_coord = ( r * cos(t) ) * sfx + cx;
	    y_coord = ( r * sin(t) ) * sfy + cy;

	    // set value of polar img
	    DIRECT_A2D_ELEM(polarImg,ri,ti) = interpolate(cartIm,x_coord,y_coord);
	}
    }


    return polarImg;
}

/* bilinear interpolation */
double ProgAngularAssignmentMag::interpolate(MultidimArray<double> &cartIm,
					     double &x_coord, double &y_coord){
    double val;
    size_t xf = floor(x_coord);
    size_t xc = ceil(x_coord);
    size_t yf = floor(y_coord);
    size_t yc = ceil(y_coord);

    if ( (xf == xc) && ( yf == yc )){
	val = dAij(cartIm, xc, yc);
    }
    else if (xf == xc){ // linear
	val = dAij(cartIm, xf, yf) + (y_coord - yf) * ( dAij(cartIm, xf, yc) - dAij(cartIm, xf, yf) );
    }
    else if(yf == yc){ // linear
	val = dAij(cartIm, xf, yf) + (x_coord - xf) * ( dAij(cartIm, xc, yf) - dAij(cartIm, xf, yf) );
    }
    else{ // bilinear
	val = ((double)(( dAij(cartIm,xf,yf)*(yc-y_coord) + dAij(cartIm,xf,yc)*(y_coord-yf) ) * (xc - x_coord)) +
	       (double)(( dAij(cartIm,xc,yf)*(yc-y_coord) + dAij(cartIm,xc,yc)*(y_coord-yf) ) * (x_coord - xf))
	       )  / (double)( (xc - xf)*(yc - yf) );
    }

    return val;

}

/* its an experiment for implement fftshift*/
void ProgAngularAssignmentMag::completeFourierShift(MultidimArray<double> &in, MultidimArray<double> &out){

    // correct output size
    out.resizeNoCopy(in);

    size_t Cf = (size_t)(YSIZE(in)/2.0 + 0.5);      //(Ydim/2.0 + 0.5);
    size_t Cc = (size_t)(XSIZE(in)/2.0 + 0.5);      //(Xdim/2.0 + 0.5);

    size_t ff, cc;
    for(size_t f = 0; f < YSIZE(in); f++){
	ff = (f + Cf) % YSIZE(in);
	for(size_t c = 0; c < XSIZE(in); c++){
	    cc = (c + Cc) % XSIZE(in);
	    DIRECT_A2D_ELEM(out, ff, cc) = DIRECT_A2D_ELEM(in,f,c);
	}
    }
}

/* its an experiment for implement fftshift*/
void ProgAngularAssignmentMag::halfFourierShift(MultidimArray<double> &in, MultidimArray<double> &out){
    size_t Cf = (size_t)(Ydim/2.0 + 0.5);
    out.resizeNoCopy(in);

    size_t ff, cc;
    for(size_t f = 0; f < Ydim; f++){
	ff = (f + Cf) % Ydim;
	for(size_t c = 0; c < Cf; c++){
	    cc = c;
	    DIRECT_A2D_ELEM(out, ff, cc) = DIRECT_A2D_ELEM(in,f,c);
	}
    }

}



/* experiment for GCC matrix product F1 .* conj(F2)
*
*/
void ProgAngularAssignmentMag::ccMatrix(MultidimArray< std::complex<double>> &F1,
                                        MultidimArray< std::complex<double>> F2,/*reference image*/
                                        MultidimArray<double> &result){


    result.resizeNoCopy(YSIZE(F1),2*(XSIZE(F1)-1));

    //    CorrelationAux aux2;
    //    correlation_matrix(F1,F2,result,aux2);

    //double mdSize=-dSize;

    CorrelationAux aux;
    aux.transformer1.setReal(result);
    aux.transformer1.setFourier(F1);
    // Multiply FFT1 .* FFT2'
    double a, b, c, d; // a+bi, c+di
    double dSize=MULTIDIM_SIZE(result);
    //    double mdSize=-dSize;

    double *ptrFFT2=(double*)MULTIDIM_ARRAY(F2);
    double *ptrFFT1=(double*)MULTIDIM_ARRAY(aux.transformer1.fFourier);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
    {
        a=(*ptrFFT1)*dSize;
        b=(*(ptrFFT1+1))*dSize;
        c=(*ptrFFT2++);
        d=(*ptrFFT2++)*(-1); //(-1);
        //        GCC
        //        *ptrFFT1++ = a*c-b*d;
        //        *ptrFFT1++ = b*c+a*d;
        // // for Compactly supported correlation
        // // F2 is reference image
        *ptrFFT1++ = (a*c-b*d)/((c*c+d*d)+0.001);
        *ptrFFT1++ = (b*c+a*d)/((c*c+d*d)+0.001);

    }

    // from correlation_matrix
    //    // Multiply FFT1 * FFT2'
    //    double mdSize=-dSize;
    //    double a, b, c, d; // a+bi, c+di
    //    double *ptrFFT2=(double*)MULTIDIM_ARRAY(FF2);
    //    double *ptrFFT1=(double*)MULTIDIM_ARRAY(FF1);
    //    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FF1)
    //    {
    //        a=*ptrFFT1++;
    //        b=*ptrFFT1++;
    //        c=(*ptrFFT2)*dSize;
    //        d=(*(ptrFFT2+1))*mdSize;
    //        *ptrFFT2++ = a*c-b*d;
    //        *ptrFFT2++ = b*c+a*d;
    //    }

    aux.transformer1.inverseFourierTransform();
    CenterFFT(result, true);
    result.setXmippOrigin();
}

/*   try PhaseCorr only for shift  ( F1 .* conj(F2) ) ./ ||  F1 .* conj(F2) ||²
*
*/
void ProgAngularAssignmentMag::ccMatrixPCO(MultidimArray< std::complex<double>> &F1,
                                        MultidimArray< std::complex<double>> &F2,
                                        MultidimArray<double> &result){


    result.resizeNoCopy(YSIZE(F1),2*(XSIZE(F1)-1));

    CorrelationAux aux;
    aux.transformer1.setReal(result);
    aux.transformer1.setFourier(F1);
    // Multiply FFT1 * FFT2'
    double a, b, c, d; // a+bi, c+di
    double *ptrFFT2=(double*)MULTIDIM_ARRAY(F2);
    double *ptrFFT1=(double*)MULTIDIM_ARRAY(aux.transformer1.fFourier);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
    {
        a=*ptrFFT1;
        b=*(ptrFFT1+1);
        c=(*ptrFFT2++);
        d=(*ptrFFT2++)*(-1);
        // phase corr only
        double den = (a*c-b*d)*(a*c-b*d) + (b*c+a*d)*(b*c+a*d);
        *ptrFFT1++ = (a*c-b*d)/(den+0.001);
        *ptrFFT1++ = (b*c+a*d)/(den+0.001);
    }

    aux.transformer1.inverseFourierTransform();
    CenterFFT(result, true);
    result.setXmippOrigin();
}


/* select n_bands of polar representation of magnitude spectrum */
void ProgAngularAssignmentMag::selectBands(MultidimArray<double> &in, MultidimArray<double> &out){

    int colStop = XSIZE(out);
    int rowStop = YSIZE(out);
    int i, j;
    // 0:179 and duplicate data
    for (i = 0; i < rowStop; i++){
	for (j = 0; j < colStop; j++){
            dAij(out,i,j) = dAij(in, startBand+i, j);
	}
    }

}

/* gets maximum value for each column*/
void ProgAngularAssignmentMag::maxByColumn(MultidimArray<double> &in,
					   MultidimArray<double> &out){

    out.resizeNoCopy(1,XSIZE(in));
    int f, c;
    double maxVal, val2;
    for(c = 0; c < XSIZE(in); c++){
	maxVal = dAij(in, 0, c);
	for(f = 1; f < YSIZE(in); f++){
	    val2 = dAij(in, f, c);
	    if (val2 > maxVal)
		maxVal = val2;
	}
	dAi(out,c) = maxVal;
    }
}

/* gets maximum value for each column*/
void ProgAngularAssignmentMag::meanByColumn(MultidimArray<double> &in,
                                           MultidimArray<double> &out){

    out.resizeNoCopy(1,XSIZE(in));
    int f, c;
    double val, val2;
    for(c = 0; c < XSIZE(in); c++){
        val = dAij(in, 0, c);
        for(f = 1; f < YSIZE(in); f++){
            val2 = dAij(in, f, c);
            val += val2;
        }
        dAi(out,c) = val/YSIZE(in);
    }
}

/* gets maximum value for each row */
void ProgAngularAssignmentMag::maxByRow(MultidimArray<double> &in,
					MultidimArray<double> &out){
    out.resizeNoCopy(1,YSIZE(in));
    int f, c;
    double maxVal, val2;
    for(f = 0; f < YSIZE(in); f++){
	maxVal = dAij(in, f, 0);
	for(c = 1; c < XSIZE(in); c++){
	    val2 = dAij(in, f, c);
	    if (val2 > maxVal)
		maxVal = val2;
	}
	dAi(out,f) = maxVal;
    }
}

/* gets maximum value for each row */
void ProgAngularAssignmentMag::meanByRow(MultidimArray<double> &in,
                                        MultidimArray<double> &out){
    out.resizeNoCopy(1,YSIZE(in));
    int f, c;
    double val, val2;
    for(f = 0; f < YSIZE(in); f++){
        val = dAij(in, f, 0);
        for(c = 1; c < XSIZE(in); c++){
            val2 = dAij(in, f, c);
            val += val2;
        }
        dAi(out,f) = val/XSIZE(in);
    }
}

/*quadratic interpolation for location of peak in crossCorr vector*/
double quadInterp(const int Idx, MultidimArray<double> &in){
    double InterpIdx = Idx - ( ( dAi(in,Idx+1) - dAi(in,Idx-1) ) / ( dAi(in,Idx+1) + dAi(in,Idx-1) - 2*dAi(in, Idx) ) )/2.;
    return InterpIdx;
}

/* Only for 180 angles */
/* approach which selects only two locations of maximum peaks in ccvRot */
void ProgAngularAssignmentMag::rotCandidates3(MultidimArray<double> &in,
					      std::vector<double> &cand,
					      const size_t &size, int *nPeaksFound){
    double max1 = -10.;
    int idx1 = 0;
    double max2 = -10.;
    int idx2 = 0;
    int i;
    int cont = 0;
    *(nPeaksFound) = cont;

    for(i = 89; i < 271; i++){ // only look for in range 90:-90
	// current value is a peak value?
	if ( (dAi(in,i) > dAi(in,i-1)) && (dAi(in,i) > dAi(in,i+1)) ){
	    cont++;
	    if( dAi(in,i) > max1){
		max2 = max1;
		idx2 = idx1;
		max1 = dAi(in,i);
		idx1 = i;
	    }
	    else if( dAi(in,i) > max2 && dAi(in,i) != max1 ){
		max2 = dAi(in,i);
		idx2 = i;
	    }
	}
    }
    int maxAccepted = 2;
    maxAccepted = ( cont < maxAccepted) ? cont : maxAccepted;
    if(cont){
	std::vector<int> temp(2,0);
        temp[0] = idx1;
        temp[1] = idx2;
        int tam = 2*maxAccepted;
	*(nPeaksFound) = tam;
	cand.reserve(tam);
	double interpIdx; // quadratic interpolated location of peak
	for(i = 0; i < maxAccepted; i++){
	    interpIdx = quadInterp(temp[i], in);
            cand[i] =  double( size - 1 )/2. - interpIdx;
	    cand[i+maxAccepted] =(cand[i]>0) ? cand[i] + 180 : cand[i] - 180 ;
	}
    }
    else{
	*(nPeaksFound) = 0;
    }
}

/* approach which selects only two locations of maximum peaks in ccvRot */
void ProgAngularAssignmentMag::rotCandidates2(MultidimArray<double> &in,
					      std::vector<double> &cand,
					      const size_t &size, int *nPeaksFound){
    const int maxNumPeaks = 20;
    double max1 = -10.;
    int idx1 = 0;
    double max2 = -10.;
    int idx2 = 0;
    int i;
    int cont = 0;
    *(nPeaksFound) = cont;
    for(i = 89/*1*/; i < 271/*size-1*/; i++){
	// current value is a peak value?
	if ( (dAi(in,i) > dAi(in,i-1)) && (dAi(in,i) > dAi(in,i+1)) ){
	    cont++;
	    if( dAi(in,i) > max1){
		max2 = max1;
		idx2 = idx1;
		max1 = dAi(in,i);
		idx1 = i;
	    }
	    else if( dAi(in,i) > max2 && dAi(in,i) != max1 ){
		max2 = dAi(in,i);
		idx2 = i;
	    }
	}
    }

    if( cont > maxNumPeaks){
	printf("reaches max number of peaks!\n");
    }

    int maxAccepted = 2;

    maxAccepted = ( cont < maxAccepted) ? cont : maxAccepted;

    if(cont){
	std::vector<int> temp(2,0);
	temp[0] = idx1;
	temp[1] = idx2;
	int tam = 2*maxAccepted;
	*(nPeaksFound) = tam;
	cand.reserve(tam);
	for(i = 0; i < maxAccepted; i++){
	    cand[i] = dAi(axRot,temp[i]);
	    cand[i+maxAccepted] =(cand[i]>0) ? cand[i] + 180 : cand[i] - 180 ;
	}

    }
    else{
	printf("no peaks found!\n");
    }
}

/* candidates to best rotation*/
void ProgAngularAssignmentMag::rotCandidates(MultidimArray<double> &in,
					     std::vector<double> &cand,
					     const size_t &size, int *nPeaksFound){
    const int maxNumPeaks = 30;
    int maxAccepted = 4;
    int *peakPos = (int*) calloc(maxNumPeaks,sizeof(int));
    int cont = 0;
    *(nPeaksFound) = cont;
    int i;
    for(i = 89/*1*/; i < 271/*size-1*/; i++){

	if ( (dAi(in,i) > dAi(in,i-1)) && (dAi(in,i) > dAi(in,i+1)) ){
	    peakPos[cont] = i;
	    cont++;
	    *(nPeaksFound) = cont;
	}
    }

    maxAccepted = ( *(nPeaksFound) < maxAccepted) ? *(nPeaksFound) : maxAccepted;

    if( *(nPeaksFound) > maxNumPeaks)
	printf("reaches max number of peaks!\n");

    if(cont){
	std::vector<int> temp(*(nPeaksFound),0);
	for(i = 0; i < *(nPeaksFound); i++){
	    temp[i] = peakPos[i];
	}
	// delete peakPos
	free(peakPos);

	// sorting first in case there are more than maxAccepted peaks
	std::sort(temp.begin(), temp.end(), [&](int i, int j){return dAi(in,i) > dAi(in,j); } );

	int tam = 2*maxAccepted; //
	*(nPeaksFound) = tam;
	cand.reserve(tam);
	for(i = 0; i < maxAccepted; i++){
	    cand[i] = dAi(axRot,temp[i]);
	    cand[i+maxAccepted] =(cand[i]>0) ? cand[i] + 180 : cand[i] - 180 ;
	}
    }
    else{
	printf("no peaks found!\n");
	// delete peakPos
	free(peakPos);
    }

}

/* instace of "delay axes" for assign rotation and traslation candidates*/
void ProgAngularAssignmentMag::_delayAxes(const size_t &Ydim, const size_t &Xdim, const size_t &n_ang){
    axRot.resize(1,1,1,n_ang);
    axTx.resize(1,1,1,Xdim);
    axTy.resize(1,1,1,Ydim);

    double M = double(n_ang - 1)/2.;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axRot){
        dAi(axRot,i) = ceil(M - i); // cambié hoy 18 junio, estaba perdiendo precisión
    }
    M = double(Xdim - 1)/2.0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axTx){
        dAi(axTx,i) = ceil(M - i);
    }
    M = double(Ydim - 1)/2.0;
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axTy){
        dAi(axTy,i) = ceil(M - i);
    }
}

/* selection of best candidate to rotation and its corresponding shift
 * called at first loop in "coarse" searching
 * shitfs are computed as maximum of CrossCorr vector
 * vector<double> cand contains candidates to relative rotation between images
*/
void ProgAngularAssignmentMag::bestCand(/*inputs*/
					MultidimArray<double> &MDaIn,
					MultidimArray< std::complex<double> > &MDaInF,
					MultidimArray<double> &MDaRef,
					std::vector<double> &cand,
					int &peaksFound,
					/*outputs*/
					double *bestCandRot,
					double *shift_x,
					double *shift_y,
					double *bestCoeff){
    *(bestCandRot) = 0;
    *(shift_x) = 0.;
    *(shift_y) = 0.;
    *(bestCoeff) = 0.0;
    double rotVar = 0.0;
    double tempCoeff;
    double tx, ty;
    MultidimArray<double> MDaRefRot;
    MultidimArray<double> MDaRefRotShift;
    MultidimArray<double> ccMatrixShift;
    MultidimArray<double> ccVectorTx;
    MultidimArray<double> ccVectorTy;
    MultidimArray< std::complex<double> > MDaRefRotF;

    MDaRefRot.setXmippOrigin();
    for(int i = 0; i < peaksFound; i++){
	rotVar = -1. * cand[i]; //
        _applyRotation(MDaRef,rotVar,MDaRefRot);
        _applyFourierImage2(MDaRefRot,MDaRefRotF);
        ccMatrix(MDaInF, MDaRefRotF, ccMatrixShift);// cross-correlation matrix / phase correlation for shift (new)
        maxByColumn(ccMatrixShift, ccVectorTx); // ccvMatrix to ccVector
        getShift(ccVectorTx,tx,XSIZE(ccMatrixShift));
        tx = -1. * tx;
        maxByRow(ccMatrixShift, ccVectorTy); // ccvMatrix to ccVector
        getShift(ccVectorTy,ty,YSIZE(ccMatrixShift));
        ty = -1. * ty;
        if ( std::abs(tx)>maxShift || std::abs(ty)>maxShift )
	    continue;
	_applyShift(MDaRefRot, tx, ty, MDaRefRotShift);
	// Pearson coeff
	pearsonCorr(MDaIn, MDaRefRotShift, tempCoeff);
	if ( tempCoeff > *(bestCoeff) ){
	    *(bestCoeff) = tempCoeff;
	    *(shift_x) = tx;
	    *(shift_y) = ty;
	    *(bestCandRot) = rotVar;
	}
    }
}

/* apply affine transform to input image avoiding some products
 *
 *                   | a b tx |
 *  affine matrix A =| c d ty |
 *                   | 0 0 1  |
 *
 */
void ProgAngularAssignmentMag::newApplyGeometry(MultidimArray<double>& __restrict__ in,
						MultidimArray<double>& __restrict__ out,
						const double &a,  const double &b,
						const double &c,  const double &d,
						const double &tx, const double &ty ){

    int nFil = YSIZE(in);
    int nCol = XSIZE(in);

    double Cx = (nCol)/2.0;
    double Cy = (nFil)/2.0;

    // constants
    double k1 = b * Cy;
    double k2 = d * Cy;

    double e1 = Cx + tx;
    double e2 = Cy + ty;

    double d1 = e1 + k1; // Cx + tx + b * Cy
    double d2 = e2 + k2; // Cy + ty + d * Cy

    double g1 = e1 - k1; // Cx + tx - b * Cy
    double g2 = e2 - k2; // Cy + ty - d * Cy


    double x1,y1,p1,q1,p2,q2,p3,q3,p4,q4;

    int x,y,rx,ry;

    int lim_x1 = 0;
    int lim_x2 = nCol-lim_x1;
    int lim_y1 = 1;
    int lim_y2 = nFil-lim_y1;

    for(x = 0.; x < Cx; x++){
	for(y = 0; y < Cy; y++){
	    x1 = a*double(x) + b*double(y);
	    y1 = c*double(x) + d*double(y);

	    // point 1 (x,y) // 4th
	    p1 = x1 + e1;
	    q1 = y1 + e2;
	    rx = x+Cx;
	    ry = y+Cy;
	    if ( (p1 > lim_x1) && (p1 < lim_x2) && (q1 > lim_y1) && (q1 < lim_y2) ){
		dAij(out, ry, rx) = interpolate( in, q1, p1);
	    }

	    // point 2 (-x, -y + Cy) // 3th
	    p2 = -x1 + d1;
	    q2 = -y1 + d2;
	    rx = -x+Cx;
	    ry = -y+2*Cy;
	    if ( (p2 > lim_x1) && (p2 < lim_x2) && (q2 > lim_y1) && (q2 < lim_y2) && (ry < lim_y2) ){
		dAij(out, ry, rx) = interpolate( in, q2, p2);
	    }

	    //point 3 (-x, -y) // 2nd
	    p3 = -x1 + e1;
	    q3 = -y1 + e2;
	    rx = -x+Cx;
	    ry = -y+Cy;
	    if ( (p3 > lim_x1) && (p3 < lim_x2) && (q3 > lim_y1) && (q3 < lim_y2) ){
		dAij(out, ry, rx) = interpolate( in, q3, p3);
	    }

	    // point 4 (x, y-Cy) // 1st
	    p4 = x1 + g1;
	    q4 = y1 + g2;
	    rx = x+Cx;
	    ry = y;
	    if ( (p4 > lim_x1) && (p4 < lim_x2) && (q4 > lim_y1) && (q4 < lim_y2) ){
		dAij(out, ry, rx) = interpolate( in, q4, p4);
	    }
	}
    }

}

/* apply rotation */
void ProgAngularAssignmentMag::_applyRotation(MultidimArray<double> &MDaRef, double &rot,
					      MultidimArray<double> &MDaRefRot){
    // Transform matrix
    Matrix2D<double> A(3,3);
    A.initIdentity();
    double ang, cosine, sine;
    ang = DEG2RAD(rot);
    cosine = cos(ang);
    sine = sin(ang);

    // rotation
    MAT_ELEM(A,0, 0) = cosine;
    MAT_ELEM(A,0, 1) = sine;
    MAT_ELEM(A,1, 0) = -sine;
    MAT_ELEM(A,1, 1) = cosine;

    // Shift
    MAT_ELEM(A,0, 2) = 0.;
    MAT_ELEM(A,1, 2) = 0.;

    applyGeometry(LINEAR, MDaRefRot, MDaRef, A, IS_NOT_INV, DONT_WRAP);

}

/* apply traslation */
void ProgAngularAssignmentMag::_applyShift(MultidimArray<double> &MDaRef,
					   double &tx, double &ty,
					   MultidimArray<double> &MDaRefShift){
    // Transform matrix
    Matrix2D<double> A(3,3);
    A.initIdentity();

    // Shift
    MAT_ELEM(A,0, 2) = tx;
    MAT_ELEM(A,1, 2) = ty;

    applyGeometry(LINEAR, MDaRefShift, MDaRef, A, IS_NOT_INV, DONT_WRAP);
}

/* finds shift as maximum of ccVector */
void ProgAngularAssignmentMag::getShift(MultidimArray<double> &ccVector, double &shift, const size_t &size){
    double maxVal = -10.;
    int idx;
    int i;
    int lb= int(size/2-maxShift);
    int hb=int(size/2+maxShift);
    for(i = lb; i < hb; i++){
	if(dAi(ccVector,i) > maxVal){
	    maxVal = dAi(ccVector,i);
	    idx = i;
	}
    }
    // interpolate value
    double interpIdx;
    interpIdx = quadInterp(idx, ccVector);
    shift = double( size - 1 )/2. - interpIdx;
}


/* finds rot as maximum of ccVector for a region near the center */
void ProgAngularAssignmentMag::getRot(MultidimArray<double> &ccVector, double &rot, const size_t &size, const double &oldAngle){
    double maxVal = -10.;
    int idx;
    int i;
    int lb= int(size/2-5);
    int hb=int(size/2+5);
    for(i = lb; i < hb+1; i++){
        if(dAi(ccVector,i) > maxVal){
            maxVal = dAi(ccVector,i);
            idx = i;
        }
    }
    // interpolate value
    double interpIdx;
    interpIdx = quadInterp(idx, ccVector);
    rot = double( size - 1 )/2. - interpIdx;
}



/* Structural similarity SSIM index Coeff */
void ProgAngularAssignmentMag::ssimIndex(MultidimArray<double> &X, MultidimArray<double> &Y, double &coeff){

    // covariance
    double X_m, Y_m, X_std, Y_std;
    double c1, c2, L;
    arithmetic_mean_and_stddev(X, X_m, X_std);
    arithmetic_mean_and_stddev(Y, Y_m, Y_std);

    double prod_mean = mean_of_products(X, Y);
    double covariace = prod_mean - (X_m * Y_m);

    L = 1;
    c1 = (0.01*L) * (0.01*L);
    c2 = (0.03*L) * (0.03*L); // estabilidad en división


    coeff = ( (2*X_m*Y_m + c1)*(2*covariace+c2) )/( (X_m*X_m + Y_m*Y_m + c1)*(X_std*X_std + Y_std*Y_std + c2) );
}

/* selection of best candidate to rotation and its corresponding shift
 * called at second loop in a little bit more strict searching
 * shitfs are computed as maximum of CrossCorr vector +0.5 / -0.5
 * vector<double> cand contains candidates to relative rotation between images with larger CrossCorr-coeff after first loop
*/
void ProgAngularAssignmentMag::bestCand2(/*inputs*/
					 MultidimArray<double> &MDaIn,
					 MultidimArray< std::complex<double> > &MDaInF,
					 MultidimArray<double> &MDaRef,
					 std::vector<double> &cand,
					 int &peaksFound,
					 /*outputs*/
					 double *bestCandRot,
					 double *shift_x,
					 double *shift_y,
					 double *bestCoeff){
    *(bestCandRot) = 0;
    *(shift_x) = 0.;
    *(shift_y) = 0.;
    *(bestCoeff) = -10.0;
    double rotVar = 0.0;
    double tempCoeff;
    double tx, ty;
    std::vector<double> vTx, vTy;
    MultidimArray<double> MDaRefRot;
    MultidimArray<double> MDaRefRotShift;
    MultidimArray<double> ccMatrixShift;
    // compare with bestCand (different initialization) line 956
    MultidimArray<double> ccVectorTx;
    MultidimArray<double> ccVectorTy;
    MultidimArray< std::complex<double> > MDaRefRotF;

    MDaRefRot.setXmippOrigin();
    for(int i = 0; i < peaksFound; i++){
	rotVar = -1. * cand[i];
	_applyRotation(MDaRef,rotVar,MDaRefRot);

	_applyFourierImage2(MDaRefRot,MDaRefRotF); // fourier --> F2_r

	ccMatrix(MDaInF, MDaRefRotF, ccMatrixShift);// cross-correlation matrix
	maxByColumn(ccMatrixShift, ccVectorTx); // ccvMatrix to ccVector
        getShift(ccVectorTx,tx,XSIZE(ccMatrixShift));
	tx = -1. * tx;
	maxByRow(ccMatrixShift, ccVectorTy); // ccvMatrix to ccVector
        getShift(ccVectorTy,ty,YSIZE(ccMatrixShift));
	ty = -1. * ty;

	if ( std::abs(tx)>maxShift || std::abs(ty)>maxShift ) // 10 es elegido pero debo poner criterio automático
	    continue;

	//*********** when strict, after first loop ***************
	// posible shifts
	vTx.push_back(tx);
	vTx.push_back(tx+0.5);
	vTx.push_back(tx-0.5);
	vTy.push_back(ty);
	vTy.push_back(ty+0.5);
	vTy.push_back(ty-0.5);

	for(int j = 0; j < 3; j++){
	    for (int k = 0; k < 3; k++){
		// translate rotated version of MDaRef
		_applyShift(MDaRefRot, vTx[j], vTy[k], MDaRefRotShift);
		// Pearson coeff

		pearsonCorr(MDaIn, MDaRefRotShift, tempCoeff);
		//        std::cout << "myCorr(f1,f2_rt): " << tempCoef << std::endl;
		if ( tempCoeff > *(bestCoeff) ){
		    *(bestCoeff) = tempCoeff;
		    *(shift_x) = vTx[j];
		    *(shift_y) = vTy[k];
		    *(bestCandRot) = rotVar;
		}
	    }
	}


    }

}

/* apply rotation */
void ProgAngularAssignmentMag::_applyRotationAndShift(MultidimArray<double> &MDaRef, double &rot, double &tx, double &ty,
						      MultidimArray<double> &MDaRefRot){
    // Transform matrix
    Matrix2D<double> A(3,3);
    A.initIdentity();
    double ang, cosine, sine;
    ang = DEG2RAD(rot);
    cosine = cos(ang);
    sine = sin(ang);

    // rotation
    MAT_ELEM(A,0, 0) = cosine;
    MAT_ELEM(A,0, 1) = sine;
    MAT_ELEM(A,1, 0) = -sine;
    MAT_ELEM(A,1, 1) = cosine;

    // Shift
    MAT_ELEM(A,0, 2) = tx;
    MAT_ELEM(A,1, 2) = ty;

    applyGeometry(LINEAR, MDaRefRot, MDaRef, A, IS_NOT_INV, DONT_WRAP);

}



