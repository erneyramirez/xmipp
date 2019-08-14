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

ProgAngularAssignmentMag::ProgAngularAssignmentMag() {
	produces_a_metadata = true;
	each_image_produces_an_output = false;
}

ProgAngularAssignmentMag::~ProgAngularAssignmentMag() {

}

void ProgAngularAssignmentMag::defineParams() {
	XmippMetadataProgram::defineParams();
	//usage
	addUsageLine(
			"Generates a list of candidates for angular assignment for each experimental image");
	//params
	//    addParamsLine("   -i <md_file>               : Metadata file with input experimental projections");
	//    addParamsLine("   -o <md_file>               : Metadata file with output projections");
	addParamsLine("   -ref <md_file>             : Metadata file with input reference projections");
	addParamsLine("  [-odir <outputDir=\".\">]   : Output directory");
	addParamsLine("  [-sym <symfile=c1>]         : Enforce symmetry in projections");
	addParamsLine("  [-sampling <sampling=1.>]         : sampling");
	addParamsLine("  [--Nsimultaneous <Nsim=1.>]       : number of simultaneous");
}

// Read arguments ==========================================================
void ProgAngularAssignmentMag::readParams() {
	XmippMetadataProgram::readParams();
	fnIn = XmippMetadataProgram::fn_in;
	fnOut = XmippMetadataProgram::fn_out;
	fnRef = getParam("-ref");
	fnDir = getParam("-odir");
	sampling = getDoubleParam("-sampling");
	XmippMetadataProgram::oroot = fnDir;
	fnSym = getParam("-sym");
	Nsim=getIntParam("--Nsimultaneous"); //unused when MPI=1
}

// Show ====================================================================
void ProgAngularAssignmentMag::show() {
	if (verbose > 0) {
		printf("%d reference images of %d x %d\n",int(sizeMdRef),int(Xdim),int(Ydim));
		printf("%d exp images of %d x %d in this group\n", int(sizeMdIn),int(Xdim), int(Ydim));
		//        printf("imgcc %d x %d from mdIn:%d, mdRef:%d\n", int(YSIZE(imgcc)), int(XSIZE(imgcc)), int(sizeMdIn), int(sizeMdRef));
		printf("\nstartBand= %d\n", int(startBand));
		printf("finalBand= %d\n", int(finalBand));
		printf("n_bands= %d\n", int(n_bands));

		XmippMetadataProgram::show();
		//        std::cout << "Input metadata              : "  << fnIn        << std::endl;
		std::cout << "Input references: " << fnRef << std::endl;
		//        std::cout << "Output directory            : "  << fnDir       << std::endl;
		//        if (fnSym != "")
		//            std::cout << "Symmetry for projections    : "  << fnSym << std::endl;
		std::cout << "sampling: " << sampling << std::endl;
	}
}

void ProgAngularAssignmentMag::startProcessing() {
	XmippMetadataProgram::startProcessing();
}

/* print in console some values of double MultidimArray */
void ProgAngularAssignmentMag::printSomeValues(MultidimArray<double> &MDa) {
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; j++)
			std::cout << "val: " << DIRECT_A2D_ELEM(MDa, i, j) << std::endl;
}

// compute variance respect to principal candidate
void neighVariance(Matrix1D<double> &neigh, double &retVal) {
	// todo check if could be some similar in my case, it seems faster
	//	mean=stddev=0;
	//    if (vdim == 0)
	//        return;
	//
	//    double sum = 0, sum2 = 0;
	//    for (size_t j = 0; j < vdim; ++j)
	//    {
	//    	double val=VEC_ELEM(*this,j);
	//        sum+=val;
	//        sum2+=val*val;
	//    }
	//    mean=sum/vdim;
	//    stddev=sum2/vdim-mean*mean;

	double val;
	double sum = 0.;
	double candVal = VEC_ELEM(neigh, 0); // value of candidate
	double diff = 0.; // difference between j and jp neighbors
	double N = neigh.vdim;
	for (int j = 0; j < N; ++j) {
		val = VEC_ELEM(neigh, j);
		diff = val - candVal;
		sum += diff * diff;
	}
	retVal = sum / N;
}

void ProgAngularAssignmentMag::computingNeighborGraph() {

	//	std::ofstream outfile("/home/jeison/Escritorio/testNeighbours.txt");
	//	outfile<< "Idx" << "\t" << "distance" <<"    \t    \n\n";

	double factor = 180. / 3.141592653;
	N_neighbors = 3; // including candidate itself
	std::vector<int> allNeighborsjp(sizeMdRef, 0); // for ordering
	std::vector<int> nearNeighbors(N_neighbors, 0);
	std::vector<double> vecNearNeighborsWeights(N_neighbors, 0.);
	Matrix1D<double> distanceToj, dirj, dirjp, nearNeighborsDist;
	printf("processing neighbors graph...\n");
	FOR_ALL_OBJECTS_IN_METADATA(mdRef)
	{
		// in case this methods works fine
		// this first lines are the same in preProcess when processing reference images
		double rotj, tiltj, psij;
		mdRef.getValue(MDL_ANGLE_ROT, rotj, __iter.objId);
		mdRef.getValue(MDL_ANGLE_TILT, tiltj, __iter.objId);
		mdRef.getValue(MDL_ANGLE_PSI, psij, __iter.objId);
		distanceToj.initZeros(sizeMdRef);
		nearNeighborsDist.initZeros(N_neighbors);
		Euler_direction(rotj, tiltj, psij, dirj);
		int jp = -1;
		for (MDIterator __iter2(mdRef); __iter2.hasNext(); __iter2.moveNext()) {
			jp += 1;
			double rotjp, tiltjp, psijp;
			mdRef.getValue(MDL_ANGLE_ROT, rotjp, __iter2.objId);
			mdRef.getValue(MDL_ANGLE_TILT, tiltjp, __iter2.objId);
			mdRef.getValue(MDL_ANGLE_PSI, psijp, __iter2.objId);

			// todo check tilt value and if it is outside some range far from tiltj
			// then dont make distance computation and put a big value, i.e. "possibly not neighbor"
			Euler_direction(rotjp, tiltjp, psijp, dirjp);
			VEC_ELEM(distanceToj,jp) = spherical_distance(dirj, dirjp);
			allNeighborsjp.at(jp) = jp;

			// todo FALTA PONER SIMETRIA
		}

		//partial sort
		// todo check all this lambda functions definition, SonarCloud gave me some hints
		std::partial_sort(allNeighborsjp.begin(),
				allNeighborsjp.begin() + N_neighbors, allNeighborsjp.end(),
				[&](int i, int j) {return VEC_ELEM(distanceToj,i) < VEC_ELEM(distanceToj,j);});
		//        // full sort
		//        std::sort(allNeighboursjp.begin(), allNeighboursjp.end(),
		//                  [&](int i, int j){return VEC_ELEM(distanceToj,i) < VEC_ELEM(distanceToj,j); });

		for (int i = 0; i < N_neighbors; ++i) {
			nearNeighbors.at(i) = allNeighborsjp.at(i); //
			VEC_ELEM(nearNeighborsDist,i) = distanceToj[nearNeighbors[i]]
					* factor; // for compute mean and std;
			//			outfile<< nearNeighbors[i] << " \t  " << VEC_ELEM(nearNeighborsDist,i) <<"  \t  ";
		}
		//		outfile<<"\n";

		//		double meanAngDist=0.;
		double varAngDist = 0.;
		//		nearNeighborsDist.computeMeanAndStddev(meanAngDist,varAngDist); //compute variance as in processImage using neighVariance method
		//		varAngDist*=varAngDist;
		neighVariance(nearNeighborsDist, varAngDist);

		for (int i = 0; i < N_neighbors; ++i) {
			vecNearNeighborsWeights.at(i) = exp(
					-0.5 * VEC_ELEM(nearNeighborsDist, i) / varAngDist);
			//			outfile<< nearNeighbors[i] << " \t  " << vecNearNeighborsWeights[i] <<"  \t  ";
		}
		//		outfile<<"\n";

		// for this __iter.objId reference image
		neighboursMatrix.push_back(nearNeighbors);
		neighboursWeights.push_back(vecNearNeighborsWeights);
	}
	//	outfile.close();
}

void ProgAngularAssignmentMag::preProcess() {
	// I want this method to run using several processors
	// it looks like runs sequentially N times when mpirun -np N is launched
	// check in mpi program. Jmendez --> just set --Nsimultaneous in highres protocol
	mdIn.read(fnIn);
	mdRef.read(fnRef);

	// size of images
	size_t Zdim, Ndim;
	getImageSize(mdIn, Xdim, Ydim, Zdim, Ndim);

	// some constants
	n_rad = size_t(Xdim / 2.);

	startBand = size_t((sampling * Xdim) / 80.); // 100
	finalBand = size_t((sampling * Xdim) / (sampling * 3)); // den MaxTargetResolution related

	n_bands = finalBand - startBand;

	n_ang = size_t(180);
	n_ang2 = 2 * n_ang;
	maxShift = .10 * Xdim; // read maxShift as input parameter    .10 * Xdim;

	// read reference images
	FileName fnImgRef;
	MDRow rowRef;
	sizeMdRef = mdRef.size();

	// how many input images
	sizeMdIn = mdIn.size();

	// reference image related
	Image<double> ImgRef;
	MultidimArray<double> MDaRef(Ydim, Xdim);
	MultidimArray<std::complex<double> > MDaRefF;
	MultidimArray<std::complex<double> > MDaRefF2;
	MultidimArray<double> MDaRefFM;
	MultidimArray<double> MDaRefFMs;
	MultidimArray<double> MDaRefFMs_polarPart(n_bands, n_ang2);
	MultidimArray<std::complex<double> > MDaRefFMs_polarF;

	size_t first = 0;
	MultidimArray<double> refPolar(n_rad, n_ang2);
	MultidimArray<std::complex<double> > MDaRefAuxF;

	//computeHann(); // UNUSED precompute Hann window
	computeCircular(); //precompute circular mask

	// for storage of rot and tilt of reference images
	referenceRot.resize(sizeMdRef);
	referenceTilt.resize(sizeMdRef);

	// try to storage all data related to reference images in memory
	printf("processing reference library...\n");
	int j = -1;
	FOR_ALL_OBJECTS_IN_METADATA(mdRef)
	{
		j += 1;
		// reading image
		mdRef.getValue(MDL_IMAGE, fnImgRef, __iter.objId);
		ImgRef.read(fnImgRef);
		MDaRef = ImgRef();
		// store to call in processImage method
		double rot, tilt, psi;
		mdRef.getValue(MDL_ANGLE_ROT, rot, __iter.objId);
		mdRef.getValue(MDL_ANGLE_TILT, tilt, __iter.objId);
		mdRef.getValue(MDL_ANGLE_PSI, psi, __iter.objId);
		referenceRot.at(j) = rot;
		referenceTilt.at(j) = tilt;
		// processing reference image
		vecMDaRef.push_back(MDaRef);
		_applyFourierImage2(MDaRef, MDaRefF);
		vecMDaRefF.push_back(MDaRefF);
		//fourier of polar image in real space
		refPolar = imToPolar(MDaRef, first, n_rad);
		// fixme add another transformer for this image size because become slower when polar images have different sizes
		// ... later this day: check if its correct because i use another transformer and it looks even slower
		_applyFourierImage2(refPolar, MDaRefAuxF, n_ang);
		vecMDaRef_polarF.push_back(MDaRefAuxF);
		// fourier of polar magnitude spectra
		transformerImage.getCompleteFourier(MDaRefF2);
		_getComplexMagnitude(MDaRefF2, MDaRefFM);
		completeFourierShift(MDaRefFM, MDaRefFMs);
		MDaRefFMs_polarPart = imToPolar(MDaRefFMs, startBand, finalBand);
		_applyFourierImage2(MDaRefFMs_polarPart, MDaRefFMs_polarF, n_ang);
		vecMDaRefFMs_polarF.push_back(MDaRefFMs_polarF);
	}

	//	// fixme could be an error to set this variables globally
	//	candidatesFirstLoop.resize(sizeMdRef);
	//	Idx.resize(sizeMdRef);
	//	candidatesFirstLoopCoeff.resize(sizeMdRef);
	//	bestTx.resize(sizeMdRef);
	//	bestTy.resize(sizeMdRef);
	//	bestPsi.resize(sizeMdRef);

	mdOut.setComment(
			"experiment for metadata output containing data for reconstruction");

	// Define the neighborhood graph
	//	computingNeighborGraph();

}

void ProgAngularAssignmentMag::processImage(const FileName &fnImg,
		const FileName &fnImgOut, const MDRow &rowIn, MDRow &rowOut) {

	// experimental image related
	rowOut = rowIn;

	// input image related
	MDRow rowRef;
	Image<double> ImgIn;
	MultidimArray<double> MDaIn(Ydim, Xdim);
	MultidimArray<std::complex<double> > MDaInF;
	MultidimArray<std::complex<double> > MDaInF2;
	MultidimArray<double> MDaInFM;
	MultidimArray<double> MDaInFMs;
	MultidimArray<double> MDaInFMs_polarPart(n_bands, n_ang2);
	MultidimArray<std::complex<double> > MDaInFMs_polarF;

	// processing input image
	ImgIn.read(fnImg);
	MDaIn = ImgIn();
	_applyFourierImage2(MDaIn, MDaInF);
	transformerImage.getCompleteFourier(MDaInF2);
	_getComplexMagnitude(MDaInF2, MDaInFM);
	completeFourierShift(MDaInFM, MDaInFMs);
	MDaInFMs_polarPart = imToPolar(MDaInFMs, startBand, finalBand);
	_applyFourierImage2(MDaInFMs_polarPart, MDaInFMs_polarF, n_ang);

	//tempCoeff = 0.0;
	double psi;
	double cc_coeff;
	double Tx;
	double Ty;
	int maxAccepted = 8;

	// todo new changing for local variables. If works, then change name
	std::vector<unsigned int> candidatesFirstLoop_local(sizeMdRef, 0);
	std::vector<unsigned int> Idx_local(sizeMdRef, 0);
	std::vector<double> candidatesFirstLoopCoeff_local(sizeMdRef, 0);
	std::vector<double> bestTx_local(sizeMdRef, 0);
	std::vector<double> bestTy_local(sizeMdRef, 0);
	std::vector<double> bestPsi_local(sizeMdRef, 0);

	MultidimArray<double> ccMatrixRot_local;
	MultidimArray<double> ccVectorRot_local;

	// loop over reference stack
	for (int k = 0; k < sizeMdRef; ++k) {
		// computing relative rotation and traslation
		ccMatrix(MDaInFMs_polarF, vecMDaRefFMs_polarF[k], ccMatrixRot_local);
		maxByColumn(ccMatrixRot_local, ccVectorRot_local);
		peaksFound = 0;
		std::vector<double> cand(maxAccepted, 0.);
		rotCandidates3(ccVectorRot_local, cand, XSIZE(ccMatrixRot_local));
		bestCand(MDaIn, MDaInF, vecMDaRef[k], cand, psi, Tx, Ty, cc_coeff);
//bestCand2(MDaIn, MDaInF, vecMDaRef[k], cand, psi, Tx, Ty,cc_coeff);
		// all results are storage for posterior partial_sort
		Idx_local[k] = k; // for sorting
		candidatesFirstLoop_local[k] = k; // for access in second loop
		candidatesFirstLoopCoeff_local[k] = cc_coeff;
		bestTx_local[k] = Tx;
		bestTy_local[k] = Ty;
		bestPsi_local[k] = psi;
	}

	/*  // skip second loop, and sort to get top candidates
	 // choose nCand of the candidates with best corrCoeff
	 int nCand = 1; // 1  3
	 std::partial_sort(Idx_local.begin(), Idx_local.begin()+nCand, Idx_local.end(),
	 [&](int i, int j){return candidatesFirstLoopCoeff_local[i] > candidatesFirstLoopCoeff_local[j]; });

	 // reading info of reference image candidate
	 double rotRef, tiltRef;
	 // reading info of reference image candidate
	 rotRef=referenceRot.at(Idx_local[0]);
	 tiltRef=referenceTilt.at(Idx_local[0]);
	 //save metadata of images with angles
	 rowOut.setValue(MDL_IMAGE,       fnImgOut);
	 rowOut.setValue(MDL_ENABLED,     1);
	 //	 rowOut.setValue(MDL_IDX,         size_t(candidatesFirstLoop[ Idx[0] ]));
	 rowOut.setValue(MDL_MAXCC,       candidatesFirstLoopCoeff_local[Idx_local[0]]);
	 rowOut.setValue(MDL_WEIGHT,      1.);
	 rowOut.setValue(MDL_WEIGHT_SIGNIFICANT,   1.);
	 rowOut.setValue(MDL_ANGLE_ROT,   rotRef);
	 rowOut.setValue(MDL_ANGLE_TILT,  tiltRef);
	 rowOut.setValue(MDL_ANGLE_PSI,   realWRAP(bestPsi_local[Idx_local[0]],-180.,180.) );
	 rowOut.setValue(MDL_SHIFT_X,     -bestTx_local[Idx_local[0]]);
	 rowOut.setValue(MDL_SHIFT_Y,     -bestTy_local[Idx_local[0]]);

	 // */

	// at this point all correlations are computed
	// now I need to check neighbors, trying to select the most probable candidate using neighbors information
	/*  // PERCENTILES
	 std::vector<double> percentVect(sizeMdRef,0);

	 std::vector<double> allVarTx(sizeMdRef);
	 std::vector<int>    IdxVarTx(sizeMdRef);
	 std::vector<double> allVarTy(sizeMdRef);
	 std::vector<int>    IdxVarTy(sizeMdRef);
	 std::vector<double> allVarPsi(sizeMdRef); // not sure if is good to check this two last variances
	 std::vector<double> allVarCC(sizeMdRef);
	 Matrix1D<double> neighTx, neighTy, neighPsi, neighCC;
	 double varTx, varTy, varPsi, varCC;
	 k=-1;
	 double dk;
	 for(int countRefImg = 0; countRefImg < sizeMdRef; countRefImg++){
	 k+=1;
	 neighTx.initZeros(N_neighbors);
	 neighTy.initZeros(N_neighbors);
	 neighPsi.initZeros(N_neighbors);
	 neighCC.initZeros(N_neighbors);
	 varTx=0.;
	 varTy=0.;
	 varPsi=0.;
	 varCC=0.;
	 double norm=0.;
	 double weight=0.;
	 int neighborIdx;
	 for(int jj=0;jj<N_neighbors;++jj){
	 neighborIdx=neighboursMatrix.at(countRefImg).at(jj);
	 weight=neighboursWeights.at(countRefImg).at(jj);
	 norm+=weight; //unused now
	 // variables values of neighbors
	 VEC_ELEM(neighTx,jj)=bestTx[neighborIdx];
	 VEC_ELEM(neighTy,jj)=bestTy[neighborIdx];
	 VEC_ELEM(neighPsi,jj)=bestRot[neighborIdx]; // if want to compute variance I have to put in correct range, i.e -90 == 270
	 VEC_ELEM(neighCC,jj)=candidatesFirstLoopCoeff[neighborIdx];

	 }

	 neighVariance(neighTx,varTx);
	 neighVariance(neighTy,varTy);
	 neighVariance(neighCC,varCC);

	 dk=k;
	 percentVect[k]=1.-dk/sizeMdRef;
	 IdxVarTx[k]=k;
	 IdxVarTy[k]=k;
	 allVarTx[k]=varTx;
	 allVarTy[k]=varTy;
	 }

	 // sort
	 std::sort(IdxVarTx.begin(), IdxVarTx.end(),
	 [&](int i, int j){return allVarTx[i] < allVarTx[j]; });
	 std::sort(IdxVarTy.begin(), IdxVarTy.end(),
	 [&](int i, int j){return allVarTy[i] < allVarTy[j]; });

	 // changing correlations
	 std::vector<double> weightedCC(sizeMdRef);
	 std::vector<int>    IdxWeightedCC(sizeMdRef);
	 k=-1;
	 for(int countRefImg = 0; countRefImg < sizeMdRef; countRefImg++){
	 k+=1;
	 double currentCC = candidatesFirstLoopCoeff[countRefImg];
	 double newCC=currentCC*percentVect[IdxVarTx[countRefImg]]*percentVect[IdxVarTy[countRefImg]];
	 IdxWeightedCC[k]=k;
	 weightedCC[k]=newCC;
	 }

	 std::sort(IdxWeightedCC.begin(), IdxWeightedCC.end(),
	 [&](int i, int j){ return weightedCC[i] > weightedCC[j]; });

	 //	for(int j=0;j<5;++j){
	 //		printf("corrBefore: %.4f \t corrAfter: %.4f \t candidate: %d \n",
	 //				candidatesFirstLoopCoeff[IdxWeightedCC[j]], weightedCC[IdxWeightedCC[j]], IdxWeightedCC[j]+1);
	 //	}
	 //	printf("\n");
	 //	std::cin.ignore();

	 double rotRef, tiltRef;
	 // reading info of reference image candidate
	 rotRef=referenceRot.at(IdxWeightedCC[0]);
	 tiltRef=referenceTilt.at(IdxWeightedCC[0]);
	 //save metadata of images with angles
	 rowOut.setValue(MDL_IMAGE,       fnImgOut);
	 rowOut.setValue(MDL_ENABLED,     1);
	 rowOut.setValue(MDL_IDX,         size_t(candidatesFirstLoop[ IdxWeightedCC[0] ]));
	 rowOut.setValue(MDL_MAXCC,       candidatesFirstLoopCoeff[IdxWeightedCC[0]]); // o weightedCC[IdxWeightedCC[0]]  candidatesFirstLoopCoeff[IdxWeightedCC[0]]
	 rowOut.setValue(MDL_WEIGHT,      1.);
	 rowOut.setValue(MDL_WEIGHT_SIGNIFICANT,   1.);
	 rowOut.setValue(MDL_ANGLE_ROT,   rotRef);
	 rowOut.setValue(MDL_ANGLE_TILT,  tiltRef);
	 rowOut.setValue(MDL_ANGLE_PSI,   bestRot[IdxWeightedCC[0]]);
	 rowOut.setValue(MDL_SHIFT_X,     -bestTx[IdxWeightedCC[0]]);
	 rowOut.setValue(MDL_SHIFT_Y,     -bestTy[IdxWeightedCC[0]]);

	 // */

	/* // SCORES
	 double fTx, fTy, fPsi, fCC, norm, weight;
	 Matrix1D<double> neighTx, neighTy, neighPsi, neighCC;
	 double varTx, varTy, varPsi, varCC, sumOfVariances;
	 k=0;
	 std::vector<unsigned int>               candidatesSecondLoop(sizeMdRef);
	 std::vector<unsigned int>               Idx2(sizeMdRef);
	 std::vector<double>                     candidatesSecondLoopScore(sizeMdRef);

	 double minSumOfVariances=1000.;
	 double newMaxCC=0.;
	 double maxScore=0.;
	 for(int countRefImg = 0; countRefImg < sizeMdRef; countRefImg++){ //
	 neighTx.initZeros(N_neighbors);
	 neighTy.initZeros(N_neighbors);
	 neighPsi.initZeros(N_neighbors);
	 neighCC.initZeros(N_neighbors);
	 varTx=0.;
	 varTy=0.;
	 varPsi=0.;
	 varCC=0.;

	 double norm=0.;
	 double weight=0.;
	 double fCC=0.;
	 int neighborIdx;

	 for(int jj=0;jj<N_neighbors;++jj){
	 neighborIdx=neighboursMatrix.at(countRefImg).at(jj);
	 weight=neighboursWeights.at(countRefImg).at(jj);
	 norm+=weight;
	 fCC+=candidatesFirstLoopCoeff[neighborIdx]*weight;
	 // todo mirar si vale la pena definir así o solo usar candidatesFirstLoopCoeff[countRefImg]
	 // todo o mirar la posibilidad de que mi vecino correle mejor porque puede que haya obtenido mejores parámetros de alineamiento y me deba quedar con esos
	 // variables values of neighbors
	 VEC_ELEM(neighTx,jj)=bestTx[neighborIdx];
	 VEC_ELEM(neighTy,jj)=bestTy[neighborIdx];
	 VEC_ELEM(neighPsi,jj)=bestRot[neighborIdx]; // if want to compute variance I have to put in correct range, i.e -90 == 270
	 VEC_ELEM(neighCC,jj)=candidatesFirstLoopCoeff[neighborIdx];
	 }

	 //score
	 fCC/=norm;
	 double score=fCC;

	 // i don't want to compute variance respect to the mean value, but respect to the value of countRefImg candidate
	 neighVariance(neighTx,varTx);
	 neighVariance(neighTy,varTy);
	 neighVariance(neighCC,varCC);

	 // looking for low variance neighborhood
	 sumOfVariances=varTx+varTy;

	 // all the results storage for sort
	 Idx2[countRefImg] = k;
	 k+=1;
	 candidatesSecondLoop[countRefImg] = countRefImg+1;
	 // new score
	 score/=sumOfVariances;

	 if((sumOfVariances<minSumOfVariances) && (fCC>newMaxCC) && (score>maxScore)){
	 candidatesSecondLoopScore[countRefImg] = score;

	 //            	printf("candidate: %d \t sOv: %.4f \t sOv_before: %.4f \t CC_now: %.4f \t CC_bef: %.4f \t score: %.4f \t fCC: %.4f\n",
	 //            			candidatesSecondLoop[countRefImg],
	 //						sumOfVariances, minSumOfVariances,
	 //						neighCC[0], newMaxCC,
	 //						score, fCC);

	 minSumOfVariances=sumOfVariances;
	 newMaxCC=fCC;
	 maxScore=score;
	 }
	 else{
	 // in this case, the score should contain the CC computed in first loop
	 candidatesSecondLoopScore[countRefImg] = candidatesFirstLoopCoeff[countRefImg];
	 }
	 }

	 // sorting based in score
	 // podria ser util definir lo de los percentiles acá también asi como en la idea "new idea"
	 // 		dk=k;
	 //      percentVect[k]=1.-dk/sizeMdRef;
	 std::sort(Idx2.begin(), Idx2.end(),
	 [&](int i, int j){ return candidatesSecondLoopScore[i] > candidatesSecondLoopScore[j]; });

	 double rotRef, tiltRef;
	 // reading info of reference image candidate
	 rotRef=referenceRot.at(Idx2[0]);
	 tiltRef=referenceTilt.at(Idx2[0]);
	 //save metadata of images with angles
	 rowOut.setValue(MDL_IMAGE,       fnImgOut);
	 rowOut.setValue(MDL_ENABLED,     1);
	 rowOut.setValue(MDL_IDX,         size_t(candidatesFirstLoop[ Idx2[0] ]));
	 rowOut.setValue(MDL_MAXCC,       candidatesFirstLoopCoeff[Idx2[0]]);
	 rowOut.setValue(MDL_WEIGHT,      1.);
	 rowOut.setValue(MDL_WEIGHT_SIGNIFICANT,   1.);
	 rowOut.setValue(MDL_ANGLE_ROT,   rotRef);
	 rowOut.setValue(MDL_ANGLE_TILT,  tiltRef);
	 rowOut.setValue(MDL_ANGLE_PSI,   bestRot[Idx2[0]]);
	 rowOut.setValue(MDL_SHIFT_X,     -bestTx[Idx2[0]]);
	 rowOut.setValue(MDL_SHIFT_Y,     -bestTy[Idx2[0]]);

	 //*/

	//	/* search rotation with polar real image representation over 10% of reference images
	int nCand = int(.10 * sizeMdRef + 1);
	std::partial_sort(Idx_local.begin(), Idx_local.begin() + nCand,
			Idx_local.end(),
			[&](int i, int j) {return candidatesFirstLoopCoeff_local[i] > candidatesFirstLoopCoeff_local[j];});

	std::vector<unsigned int> candidatesSecondLoop(nCand, 0);
	std::vector<unsigned int> Idx2(nCand, 0);
	std::vector<double> candidatesSecondLoopCoeff(nCand, 0.);
	std::vector<double> bestTx2(nCand, 0.);
	std::vector<double> bestTy2(nCand, 0.);
	std::vector<double> bestPsi2(nCand, 0.);

	// variables new approach
	size_t first = 0;
	MultidimArray<double> inPolar(n_rad, n_ang2);

	MultidimArray<double> MDaExpShiftRot2; // transform experimental
	MDaExpShiftRot2.setXmippOrigin();

	MultidimArray<double> ccMatrixRot2;
	MultidimArray<double> ccVectorRot2;

	MultidimArray<std::complex<double> > MDaInAuxF;

	for (int k = 0; k < nCand; ++k) {
		//apply transform to experimental image
		double rotVal = -1. * bestPsi_local[Idx_local[k]]; // -1. because i return parameters for reference image instead of experimental
		double trasXval = -1. * bestTx_local[Idx_local[k]];
		double trasYval = -1. * bestTy_local[Idx_local[k]];
		_applyShiftAndRotation(MDaIn, rotVal, trasXval, trasYval,
				MDaExpShiftRot2);
		//fourier experimental image
		_applyFourierImage2(MDaExpShiftRot2, MDaInF);
		// polar experimental image
		inPolar = imToPolar(MDaExpShiftRot2, first, n_rad); // this process can be expensive
		_applyFourierImage2(inPolar, MDaInAuxF, n_ang); // todo check!! , at least first time, i am using the same transformer that i use in the past loop which have a different size

		// find rotation and shift
		ccMatrix(MDaInAuxF,vecMDaRef_polarF[candidatesFirstLoop_local[Idx_local[k]]],ccMatrixRot2);
		maxByColumn(ccMatrixRot2, ccVectorRot2);
		peaksFound = 0;
		std::vector<double> cand(maxAccepted, 0.);
		rotCandidates3(ccVectorRot2, cand, XSIZE(ccMatrixRot2));
		bestCand(MDaExpShiftRot2, MDaInF, vecMDaRef[Idx_local[k]], cand, psi,Tx, Ty, cc_coeff); //
//bestCand2(MDaExpShiftRot2, MDaInF, vecMDaRef[Idx_local[k]],cand, psi, Tx, Ty, cc_coeff); //

		// if its better and shifts are within then update
		double testShiftTx = bestTx_local[Idx_local[k]] + Tx;
		double testShiftTy = bestTy_local[Idx_local[k]] + Ty;
		if (cc_coeff >= candidatesFirstLoopCoeff_local[Idx_local[k]]
				&& std::abs(testShiftTx) < maxShift
				&& std::abs(testShiftTy) < maxShift && std::abs(psi) < 2.) { // todo must define which value and why 2. 5. ?
			Idx2[k] = k;
			candidatesSecondLoop[k] = candidatesFirstLoop_local[Idx_local[k]];
			candidatesSecondLoopCoeff[k] = cc_coeff;
			bestTx2[k] = testShiftTx;
			bestTy2[k] = testShiftTy;
			bestPsi2[k] = bestPsi_local[Idx_local[k]] + psi;
		} else { //COSS suggests to degrade this candidate (i.e. coeff--> 0), because it doesn't improve
			Idx2[k] = k;
			candidatesSecondLoop[k] = candidatesFirstLoop_local[Idx_local[k]];
			candidatesSecondLoopCoeff[k] =
			candidatesFirstLoopCoeff_local[Idx_local[k]]; //candidatesFirstLoopCoeff[Idx[k]];
			// todo check if works better. It seems that for several images in VIRUS there is no any improvement, i.e. all coeff=0
			bestTx2[k] = bestTx_local[Idx_local[k]];
			bestTy2[k] = bestTy_local[Idx_local[k]];
			bestPsi2[k] = bestPsi_local[Idx_local[k]];
		}
	}

	// choose nCand of the candidates with best corrCoeff
	//	int nCand = 1; // 1  3
	std::sort(Idx2.begin(), Idx2.end(),
			[&](int i, int j) {return candidatesSecondLoopCoeff[i] > candidatesSecondLoopCoeff[j];});

	// reading info of reference image candidate
	double rotRef, tiltRef;
	// reading info of reference image candidate
	rotRef = referenceRot.at(candidatesSecondLoop[Idx2[0]]);
	tiltRef = referenceTilt.at(candidatesSecondLoop[Idx2[0]]);
	//save metadata of images with angles
	rowOut.setValue(MDL_IMAGE, fnImgOut);
	rowOut.setValue(MDL_ENABLED, 1);
	//rowOut.setValue(MDL_IDX,         size_t(candidatesSecondLoop[ Idx2[0] ]));
	rowOut.setValue(MDL_MAXCC, candidatesSecondLoopCoeff[Idx2[0]]);
	rowOut.setValue(MDL_WEIGHT, 1.);
	rowOut.setValue(MDL_WEIGHT_SIGNIFICANT, 1.);
	rowOut.setValue(MDL_ANGLE_ROT, rotRef);
	rowOut.setValue(MDL_ANGLE_TILT, tiltRef);
	rowOut.setValue(MDL_ANGLE_PSI, realWRAP(bestPsi2[Idx2[0]], -180., 180.));
	rowOut.setValue(MDL_SHIFT_X, -bestTx2[Idx2[0]]);
	rowOut.setValue(MDL_SHIFT_Y, -bestTy2[Idx2[0]]);

	// */

	//    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	//    std::cout << "Operation took "<< duration*1000 << "milliseconds" << std::endl;
}

void ProgAngularAssignmentMag::postProcess() {

	// from angularContinousAssign2
	MetaData &ptrMdOut = *getOutputMd();
	ptrMdOut.removeDisabled();
	/*  // todo ask if this is the correct way to do this weighting, because later in highres protocol there is another weighting step
	 double maxCC = -1.;
	 FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
	 {
	 double thisMaxCC;
	 ptrMdOut.getValue(MDL_MAXCC, thisMaxCC, __iter.objId);
	 if (thisMaxCC > maxCC)
	 maxCC = thisMaxCC;
	 if (thisMaxCC == 0)
	 ptrMdOut.removeObject(__iter.objId);
	 }
	 FOR_ALL_OBJECTS_IN_METADATA(ptrMdOut)
	 {
	 double thisMaxCC;
	 ptrMdOut.getValue(MDL_MAXCC, thisMaxCC, __iter.objId);
	 ptrMdOut.setValue(MDL_WEIGHT, thisMaxCC / maxCC, __iter.objId);
	 ptrMdOut.setValue(MDL_WEIGHT_SIGNIFICANT, thisMaxCC / maxCC,
	 __iter.objId);
	 } // */

	ptrMdOut.write(XmippMetadataProgram::fn_out.replaceExtension("xmd"));
}

/* Pearson Coeff. ZNCC zero-mean normalized cross-corr*/
void ProgAngularAssignmentMag::pearsonCorr(const MultidimArray<double> &X,
		MultidimArray<double> &Y, double &coeff) {

	// covariance
	double X_m, Y_m, X_std, Y_std;
	arithmetic_mean_and_stddev(X, X_m, X_std);
	arithmetic_mean_and_stddev(Y, Y_m, Y_std);

	double mean_prod = mean_of_products(X, Y);
	double covariace = mean_prod - (X_m * Y_m);

	coeff = covariace / (X_std * Y_std);
}

/* Arithmetic mean and stdDev for Pearson Coeff */
void ProgAngularAssignmentMag::arithmetic_mean_and_stddev(
		const MultidimArray<double> &data, double &avg, double &stddev) {
	data.computeAvgStdev(avg, stddev);
}

/* Arithmetic mean and stdDev for Pearson Coeff */
void ProgAngularAssignmentMag::arithmetic_mean_and_stddev(
		MultidimArray<double> &data, double &avg, double &stddev) {
	data.computeAvgStdev(avg, stddev);
}

/* Mean of products for Pearson Coeff */
double ProgAngularAssignmentMag::mean_of_products(
		const MultidimArray<double> &data1, MultidimArray<double> &data2) {
	double total = 0;
	for (int f = 0; f < Ydim; f++) {
		for (int c = 0; c < Xdim; c++) {
			total += DIRECT_A2D_ELEM(data1,f,c) * DIRECT_A2D_ELEM(data2, f, c);
		}
	}
	return total / (Xdim * Ydim);
}

/* Mean of products for Pearson Coeff */
double ProgAngularAssignmentMag::mean_of_products(MultidimArray<double> &data1,
		MultidimArray<double> &data2) {
	double total = 0;
	for (int f = 0; f < Ydim; f++) {
		for (int c = 0; c < Xdim; c++) {
			total += DIRECT_A2D_ELEM(data1,f,c) * DIRECT_A2D_ELEM(data2, f, c);
		}
	}
	return total / (Xdim * Ydim);
}

/* Normalized cross correlation*/
void ::ProgAngularAssignmentMag::normalized_cc(MultidimArray<double> &X,
		MultidimArray<double> &Y, double &value) {
	double prodXY = 0;
	double prodXX = 0;
	double prodYY = 0;
	MultidimArray<double> X2, Y2;
	X2.resizeNoCopy(X);
	Y2.resizeNoCopy(Y);

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(X)
	{
		prodXY += dAij(X,i,j) * dAij(Y, i, j);
		prodXX += dAij(X,i,j) * dAij(X, i, j);
		prodYY += dAij(Y,i,j) * dAij(Y, i, j);
	}
	double den = prodXX * prodYY;
	if (den <= 0)
		std::cout << "zero/negative denominator!!\n" << std::endl;
	else
		value = prodXY / sqrt(den);
}

/* Normalized cross correlation*/
void ::ProgAngularAssignmentMag::normalized_cc(const MultidimArray<double> &X,
		MultidimArray<double> &Y, double &value) {
	double prodXY = 0;
	double prodXX = 0;
	double prodYY = 0;
	MultidimArray<double> X2, Y2;
	X2.resizeNoCopy(X);
	Y2.resizeNoCopy(Y);

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(X)
	{
		prodXY += dAij(X,i,j) * dAij(Y, i, j);
		prodXX += dAij(X,i,j) * dAij(X, i, j);
		prodYY += dAij(Y,i,j) * dAij(Y, i, j);
	}
	double den = prodXX * prodYY;
	if (den <= 0)
		std::cout << "zero/negative denominator!!\n" << std::endl;
	else
		value = prodXY / sqrt(den);
}

/* mixed between IMED and NCC --> IMNCC ----------------------------------------------- */
void ProgAngularAssignmentMag::imNormalized_cc(const MultidimArray<double>& I1,
		const MultidimArray<double>& I2, double &value) {
	// [x,y]=meshgrid([-3:1:3],[-3:1:3])
	// format long	// w=1/sqrt(2*pi)*exp(-0.5*(x.*x+y.*y))
	double *refW;
	double w[49] = { 0.000049233388666, 0.000599785460091, 0.002688051941039,
			0.004431848411938, 0.002688051941039, 0.000599785460091,
			0.000049233388666, 0.000599785460091, 0.007306882745281,
			0.032747176537767, 0.053990966513188, 0.032747176537767,
			0.007306882745281, 0.000599785460091, 0.002688051941039,
			0.032747176537767, 0.146762663173740, 0.241970724519143,
			0.146762663173740, 0.032747176537767, 0.002688051941039,
			0.004431848411938, 0.053990966513188, 0.241970724519143,
			0.398942280401433, 0.241970724519143, 0.053990966513188,
			0.004431848411938, 0.002688051941039, 0.032747176537767,
			0.146762663173740, 0.241970724519143, 0.146762663173740,
			0.032747176537767, 0.002688051941039, 0.000599785460091,
			0.007306882745281, 0.032747176537767, 0.053990966513188,
			0.032747176537767, 0.007306882745281, 0.000599785460091,
			0.000049233388666, 0.000599785460091, 0.002688051941039,
			0.004431848411938, 0.002688051941039, 0.000599785460091,
			0.000049233388666 };

	int imiddle = YSIZE(I1) / 2;
	int jmiddle = XSIZE(I1) / 2;
	int R2max = imiddle * imiddle;
	int ysize = (int) YSIZE(I1);
	int xsize = (int) XSIZE(I1);

	double numXY = 0.;
	double denXX = 0.;
	double denYY = 0.;

	MultidimArray<double> prodImageXY = I1 * I2;
	MultidimArray<double> prodImageXX = I1 * I1;
	MultidimArray<double> prodImageYY = I2 * I2;

	for (int i = 3; i < ysize - 3; ++i) {
		int i2 = (i - imiddle) * (i - imiddle);
		for (int j = 3; j < xsize - 3; ++j) {
			int j2 = (j - jmiddle) * (j - jmiddle);
			if (i2 + j2 > R2max) // Measure only within the maximum circle
				continue;

			double prodNumXYi = DIRECT_A2D_ELEM(prodImageXY, i, j);
			double prodDenXXi = DIRECT_A2D_ELEM(prodImageXX, i, j);
			double prodDenYYi = DIRECT_A2D_ELEM(prodImageYY, i, j);
			int index = 0;
			for (int ii = -3; ii <= 3; ++ii) {
				refW = &w[index];
				index = index + 7;
				// numerator
				double *prodNumXYj = &DIRECT_A2D_ELEM(prodImageXY, i + ii,
						j - 3);
				double prodNumXYAux = (*refW) * (*prodNumXYj++);
				// Denominator XX
				double *prodDenXXj = &DIRECT_A2D_ELEM(prodImageXX, i + ii,
						j - 3);
				double prodDenXXAux = (*refW) * (*prodDenXXj++);
				// Denominator YY
				double *prodDenYYj = &DIRECT_A2D_ELEM(prodImageYY, i + ii,
						j - 3);
				double prodDenYYAux = (*refW++) * (*prodDenYYj++); //increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				numXY += prodNumXYAux * prodNumXYi;
				denXX += prodDenXXAux * prodDenXXi;
				denYY += prodDenYYAux * prodDenYYi;
			}
		}
	}

	double denom = denXX * denYY;
	value = 0;
	if (denom <= 0)
		std::cout << "zero/negative denominator!" << std::endl;
	else
		value = numXY / sqrt(denom);

}

/* mixed between IMED and ZNCC --> IMNCC ----------------------------------------------- */
void ProgAngularAssignmentMag::imZNCC(const MultidimArray<double>& I1,
		const MultidimArray<double>& I2, double &value) {

	// [x,y]=meshgrid([-3:1:3],[-3:1:3])
	// format long	// w=1/sqrt(2*pi)*exp(-0.5*(x.*x+y.*y))
	double *refW;
	double w[49] = { 0.000049233388666, 0.000599785460091, 0.002688051941039,
			0.004431848411938, 0.002688051941039, 0.000599785460091,
			0.000049233388666, 0.000599785460091, 0.007306882745281,
			0.032747176537767, 0.053990966513188, 0.032747176537767,
			0.007306882745281, 0.000599785460091, 0.002688051941039,
			0.032747176537767, 0.146762663173740, 0.241970724519143,
			0.146762663173740, 0.032747176537767, 0.002688051941039,
			0.004431848411938, 0.053990966513188, 0.241970724519143,
			0.398942280401433, 0.241970724519143, 0.053990966513188,
			0.004431848411938, 0.002688051941039, 0.032747176537767,
			0.146762663173740, 0.241970724519143, 0.146762663173740,
			0.032747176537767, 0.002688051941039, 0.000599785460091,
			0.007306882745281, 0.032747176537767, 0.053990966513188,
			0.032747176537767, 0.007306882745281, 0.000599785460091,
			0.000049233388666, 0.000599785460091, 0.002688051941039,
			0.004431848411938, 0.002688051941039, 0.000599785460091,
			0.000049233388666 };

	int imiddle = YSIZE(I1) / 2;
	int jmiddle = XSIZE(I1) / 2;
	int R2max = imiddle * imiddle;
	int ysize = (int) YSIZE(I1);
	int xsize = (int) XSIZE(I1);

	double numXY = 0.;
	double denXX = 0.;
	double denYY = 0.;

	//compute mean
	double avgI1 = I1.computeAvg();
	double avgI2 = I2.computeAvg();

	MultidimArray<double> centeredI1 = I1 - avgI1;
	MultidimArray<double> centeredI2 = I2 - avgI2;

	MultidimArray<double> prodImageXY = centeredI1 * centeredI2;
	MultidimArray<double> prodImageXX = centeredI1 * centeredI1;
	MultidimArray<double> prodImageYY = centeredI2 * centeredI2;

	for (int i = 3; i < ysize - 3; ++i) {
		int i2 = (i - imiddle) * (i - imiddle);
		for (int j = 3; j < xsize - 3; ++j) {
			int j2 = (j - jmiddle) * (j - jmiddle);
			if (i2 + j2 > R2max) // Measure only within the maximum circle
				continue;
			//        	 /* using one loop
			double prodNumXYi = DIRECT_A2D_ELEM(prodImageXY, i, j);
			double prodDenXXi = DIRECT_A2D_ELEM(prodImageXX, i, j);
			double prodDenYYi = DIRECT_A2D_ELEM(prodImageYY, i, j);
			int index = 0;
			for (int ii = -3; ii <= 3; ++ii) {
				refW = &w[index];
				index = index + 7;
				// numerator
				double *prodNumXYj = &DIRECT_A2D_ELEM(prodImageXY, i + ii,
						j - 3);
				double prodNumXYAux = (*refW) * (*prodNumXYj++);
				// Denominator XX
				double *prodDenXXj = &DIRECT_A2D_ELEM(prodImageXX, i + ii,
						j - 3);
				double prodDenXXAux = (*refW) * (*prodDenXXj++);
				// Denominator YY
				double *prodDenYYj = &DIRECT_A2D_ELEM(prodImageYY, i + ii,
						j - 3);
				double prodDenYYAux = (*refW++) * (*prodDenYYj++); //increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				prodNumXYAux += (*refW) * (*prodNumXYj++);
				prodDenXXAux += (*refW) * (*prodDenXXj++);
				prodDenYYAux += (*refW++) * (*prodDenYYj++); // increment

				numXY += prodNumXYAux * prodNumXYi;
				denXX += prodDenXXAux * prodDenXXi;
				denYY += prodDenYYAux * prodDenYYi;
			} // */
		}
	}

	double denom = denXX * denYY;
	value = 0;
	if (denom <= 0)
		std::cout << "zero/negative denominator!" << std::endl;
	else
		value = numXY / sqrt(denom);
}

void ProgAngularAssignmentMag::_applyCircularMask(
		const MultidimArray<double> &in, MultidimArray<double> &out) {

	double Cf = (Ydim + (Ydim % 2)) / 2.0; // for odd/even images
	double Cc = (Xdim + (Xdim % 2)) / 2.0;
	int pixReduc = 2;
	double rad2 = (Cf - pixReduc) * (Cf - pixReduc);
	double val = 0;
	out.initZeros(Ydim, Xdim);
	for (size_t f = 0; f < Ydim; f++) {
		for (size_t c = 0; c < Xdim; c++) {
			val = (f - Cf) * (f - Cf) + (c - Cc) * (c - Cc);
			if (val < rad2)
				DIRECT_A2D_ELEM(out, f, c) = DIRECT_A2D_ELEM(in, f, c);
		}
	}
}

/* writing out some data to file with an specified size*/
void ProgAngularAssignmentMag::_writeTestFile(MultidimArray<double> &data,
		const char* fileName, size_t nFil, size_t nCol) {
	std::ofstream outFile(fileName);
	for (int f = 0; f < nFil; f++) {
		for (int c = 0; c < nCol; c++) {
			outFile << DIRECT_A2D_ELEM(data, f, c) << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
}

/* writing out some data to file with an specified size*/
void ProgAngularAssignmentMag::_writeTestFile(const MultidimArray<double> &data,
		const char* fileName, size_t nFil, size_t nCol) {
	std::ofstream outFile(fileName);
	for (int f = 0; f < nFil; f++) {
		for (int c = 0; c < nCol; c++) {
			outFile << DIRECT_A2D_ELEM(data, f, c) << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
}

/* writing out some data to file Ydim x Xdim size*/
void ProgAngularAssignmentMag::_writeTestFile(MultidimArray<double> &data,
		const char* fileName) {
	std::ofstream outFile(fileName);
	for (int f = 0; f < Ydim; f++) {
		for (int c = 0; c < Xdim; c++) {
			outFile << DIRECT_A2D_ELEM(data, f, c) << "\t";
		}
		outFile << "\n";
	}
	outFile.close();
}

/* get COMPLETE fourier spectrum of Images. It should be changed for half */
void ProgAngularAssignmentMag::_applyFourierImage(MultidimArray<double> &data,
		MultidimArray<std::complex<double> > &FourierData) {
	transformerImage.completeFourierTransform(data, FourierData);
}

/* get COMPLETE fourier spectrum of polarRepresentation of Magnitude. It should be changed for half */
void ProgAngularAssignmentMag::_applyFourierImage(MultidimArray<double> &data,
		MultidimArray<std::complex<double> > &FourierData, const size_t &ang) {
	transformerPolarImage.completeFourierTransform(data, FourierData);
}

/*first try in using only one half of Fourier space*/
void ProgAngularAssignmentMag::_applyFourierImage2(MultidimArray<double> &data,
		MultidimArray<std::complex<double> > &FourierData) {
	transformerImage.FourierTransform(data, FourierData, true);
}

/* first try one half of fourier spectrum of polarRepresentation of Magnitude*/
void ProgAngularAssignmentMag::_applyFourierImage2(MultidimArray<double> &data,
		MultidimArray<std::complex<double> > &FourierData, const size_t &ang) {
	transformerPolarImage.FourierTransform(data, FourierData, true); // false --> true para generar copia
}

/* first try one half of fourier spectrum of polarRepresentation of image in real space*/
void ProgAngularAssignmentMag::_applyFourierImage3(MultidimArray<double> &data,
		MultidimArray<std::complex<double> > &FourierData, const size_t &ang) {
	transformerPolarRealSpace.FourierTransform(data, FourierData, true); //
}

/* get magnitude of fourier spectrum */
void ProgAngularAssignmentMag::_getComplexMagnitude(
		MultidimArray<std::complex<double> > &FourierData,
		MultidimArray<double> &FourierMag) {
	FFT_magnitude(FourierData, FourierMag);
}

/* cartImg contains cartessian  grid representation of image,
 *  rad and ang are the number of radius and angular elements*/
MultidimArray<double> ProgAngularAssignmentMag::imToPolar(
		MultidimArray<double> &cartIm, size_t &start, size_t &end) {

	int thisNbands = end - start;
	MultidimArray<double> polarImg(thisNbands, n_ang2);
	float pi = 3.141592653;
	// coordinates of center
	//    double cy = (Ydim+1)/2.0;
	//    double cx = (Xdim+1)/2.0;
	double cy = (Ydim + (Ydim % 2)) / 2.0; // for odd/even images
	double cx = (Xdim + (Xdim % 2)) / 2.0;

	// scale factors
	double sfy = (Ydim - 1) / 2.0;
	double sfx = (Xdim - 1) / 2.0;

	double delR = (double) (1.0 / (n_rad)); // n_rad-1
	double delT = 2.0 * pi / n_ang2;

	// loop through rad and ang coordinates
	double r, t, x_coord, y_coord;
	for (size_t ri = start; ri < end; ++ri) {
		for (size_t ti = 0; ti < n_ang2; ++ti) {
			r = ri * delR;
			t = ti * delT;
			x_coord = (r * cos(t)) * sfx + cx;
			y_coord = (r * sin(t)) * sfy + cy;
			// set value of polar img
			DIRECT_A2D_ELEM(polarImg,ri-start,ti) = interpolate(cartIm, x_coord,
					y_coord);
		}
	}

	//    printf("termina polar\n r, t, xcoord, ycoord = %.2f, %.2f, %.2f, %.2f\n", r, t, x_coord, y_coord);

	return polarImg;
}

/* cartImg contains cartessian  grid representation of image,
 *  rad and ang are the number of radius and angular elements
 *  this function was built for half representation of Fourier spectrum*/
MultidimArray<double> ProgAngularAssignmentMag::imToPolar2(
		MultidimArray<double> &cartIm, const size_t &rad, const size_t &ang) {
	MultidimArray<double> polarImg(rad, ang);
	float pi = 3.141592653;
	// coordinates of center
	double cy = 0.5; //(Ydim+1)/2.0;
	double cx = (Xdim + 1) / 2.0;
	// scale factors
	double sfy = (Ydim - 1) / 2.0;
	double sfx = (Xdim - 1) / 2.0;

	double delR = (double) (1.0 / (rad - 1));
	double delT = pi / ang;

	// loop through rad and ang coordinates
	double r, t, x_coord, y_coord;
	for (size_t ri = 0; ri < rad; ++ri) {
		for (size_t ti = 0; ti < ang; ++ti) {
			r = ri * delR;
			t = ti * delT;
			x_coord = (r * cos(t)) * sfx + cx;
			y_coord = (r * sin(t)) * sfy + cy;

			// set value of polar img
			DIRECT_A2D_ELEM(polarImg,ri,ti) = interpolate(cartIm, x_coord,
					y_coord);
		}
	}

	return polarImg;
}

/* bilinear interpolation */
double ProgAngularAssignmentMag::interpolate(MultidimArray<double> &cartIm,
		double &x_coord, double &y_coord) {
	double val;
	size_t xf = floor(x_coord);
	size_t xc = ceil(x_coord);
	size_t yf = floor(y_coord);
	size_t yc = ceil(y_coord);

	if ((xf == xc) && (yf == yc)) {
		val = dAij(cartIm, xc, yc);
	} else if (xf == xc) { // linear
		val = dAij(cartIm, xf, yf)
				+ (y_coord - yf)
						* ( dAij(cartIm, xf, yc) - dAij(cartIm, xf, yf));
	} else if (yf == yc) { // linear
		val = dAij(cartIm, xf, yf)
				+ (x_coord - xf)
						* ( dAij(cartIm, xc, yf) - dAij(cartIm, xf, yf));
	} else { // bilinear
		val =
				((double) (( dAij(cartIm,xf,yf) * (yc - y_coord)
						+ dAij(cartIm,xf,yc) * (y_coord - yf)) * (xc - x_coord))
						+ (double) (( dAij(cartIm,xc,yf) * (yc - y_coord)
								+ dAij(cartIm,xc,yc) * (y_coord - yf))
								* (x_coord - xf)))
						/ (double) ((xc - xf) * (yc - yf));
	}

	return val;

}

/* its an experiment for implement fftshift*/
void ProgAngularAssignmentMag::completeFourierShift(MultidimArray<double> &in,
		MultidimArray<double> &out) {

	// correct output size
	out.resizeNoCopy(in);

	size_t Cf = (size_t) (YSIZE(in) / 2.0 + 0.5); //(Ydim/2.0 + 0.5); //FIXME I think this does not work properly for even/odd images
	size_t Cc = (size_t) (XSIZE(in) / 2.0 + 0.5); //(Xdim/2.0 + 0.5); // (size_t)( 120 / 2.0 + 0.5) == (size_t)( 119 / 2.0 + 0.5)

	// TODO check this and see how to set proper value for size_t ff, cc
	//	double Cf = (YSIZE(in) + (YSIZE(in) % 2)) / 2.0; // for odd/even images
	//	double Cc = (XSIZE(in) + (XSIZE(in) % 2)) / 2.0;

	// TODO check change coordinates as (-1)^k

	size_t ff, cc;
	for (size_t f = 0; f < YSIZE(in); f++) {
		ff = (f + Cf) % YSIZE(in);
		for (size_t c = 0; c < XSIZE(in); c++) {
			cc = (c + Cc) % XSIZE(in);
			DIRECT_A2D_ELEM(out, ff, cc) = DIRECT_A2D_ELEM(in, f, c);
		}
	}
}

/* its an experiment for implement fftshift*/
void ProgAngularAssignmentMag::halfFourierShift(MultidimArray<double> &in,
		MultidimArray<double> &out) {
	size_t Cf = (size_t) (Ydim / 2.0 + 0.5);
	out.resizeNoCopy(in);

	size_t ff, cc;
	for (size_t f = 0; f < Ydim; f++) {
		ff = (f + Cf) % Ydim;
		for (size_t c = 0; c < Cf; c++) {
			cc = c;
			DIRECT_A2D_ELEM(out, ff, cc) = DIRECT_A2D_ELEM(in, f, c);
		}
	}

}

/* experiment for GCC matrix product F1 .* conj(F2)
 *
 */
void ProgAngularAssignmentMag::ccMatrix(
		const MultidimArray<std::complex<double>> &F1,
		const MultidimArray<std::complex<double>> &F2,/*reference image*/
		MultidimArray<double> &result) {

	result.resizeNoCopy(YSIZE(F1), 2 * (XSIZE(F1) - 1));

	//    CorrelationAux aux2;
	//	correlation_matrix(F1,F2,result,aux2);

	//double mdSize=-dSize;

	CorrelationAux aux;
	aux.transformer1.setReal(result);
	aux.transformer1.setFourier(F1);
	// Multiply FFT1 .* FFT2'
	double a, b, c, d; // a+bi, c+di
	double dSize = MULTIDIM_SIZE(result);
	//    double mdSize=-dSize;

	double *ptrFFT2 = (double*) MULTIDIM_ARRAY(F2);
	double *ptrFFT1 = (double*) MULTIDIM_ARRAY(aux.transformer1.fFourier);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
	{
		a = (*ptrFFT1) * dSize;
		b = (*(ptrFFT1 + 1)) * dSize;
		c = (*ptrFFT2++);
		d = (*ptrFFT2++) * (-1); //(-1);
		//        //GCC
		//        *ptrFFT1++ = a*c-b*d;
		//        *ptrFFT1++ = b*c+a*d;
		//Compactly supported correlation. F2 is reference image
		*ptrFFT1++ = (a * c - b * d) / ((c * c + d * d) + 0.001);
		*ptrFFT1++ = (b * c + a * d) / ((c * c + d * d) + 0.001);
		//        //phase correlation only
		//        double den = (a*c-b*d)*(a*c-b*d) + (b*c+a*d)*(b*c+a*d);
		//        *ptrFFT1++ = (a*c-b*d)/(den+0.001);
		//        *ptrFFT1++ = (b*c+a*d)/(den+0.001);

	}
	aux.transformer1.inverseFourierTransform();
	CenterFFT(result, true);
	result.setXmippOrigin();
}

/*   try PhaseCorr only for shift  ( F1 .* conj(F2) ) ./ ||  F1 .* conj(F2) ||²
 *
 */
void ProgAngularAssignmentMag::ccMatrixPCO(
		MultidimArray<std::complex<double>> &F1,
		MultidimArray<std::complex<double>> &F2,
		MultidimArray<double> &result) {

	result.resizeNoCopy(YSIZE(F1), 2 * (XSIZE(F1) - 1));

	CorrelationAux aux;
	aux.transformer1.setReal(result);
	aux.transformer1.setFourier(F1);
	// Multiply FFT1 * FFT2'
	double a, b, c, d; // a+bi, c+di
	double *ptrFFT2 = (double*) MULTIDIM_ARRAY(F2);
	double *ptrFFT1 = (double*) MULTIDIM_ARRAY(aux.transformer1.fFourier);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
	{
		a = *ptrFFT1;
		b = *(ptrFFT1 + 1);
		c = (*ptrFFT2++);
		d = (*ptrFFT2++) * (-1);
		// phase corr only
		double den = (a * c - b * d) * (a * c - b * d)
				+ (b * c + a * d) * (b * c + a * d);
		*ptrFFT1++ = (a * c - b * d) / (den + 0.001);
		*ptrFFT1++ = (b * c + a * d) / (den + 0.001);
	}

	aux.transformer1.inverseFourierTransform();
	CenterFFT(result, true);
	result.setXmippOrigin();
}

/* select n_bands of polar representation of magnitude spectrum */
void ProgAngularAssignmentMag::selectBands(MultidimArray<double> &in,
		MultidimArray<double> &out) {

	int colStop = XSIZE(out);
	int rowStop = YSIZE(out);
	int i, j;
	// 0:179 and duplicate data
	for (i = 0; i < rowStop; ++i) {
		for (j = 0; j < colStop; j++) {
			dAij(out,i,j) = dAij(in, startBand + i, j);
		}
	}

}

/* gets maximum value for each column*/
void ProgAngularAssignmentMag::maxByColumn(MultidimArray<double> &in,
		MultidimArray<double> &out) {

	out.resizeNoCopy(1, XSIZE(in));
	int f, c;
	double maxVal, val2;
	for (c = 0; c < XSIZE(in); c++) {
		maxVal = dAij(in, 0, c);
		for (f = 1; f < YSIZE(in); f++) {
			val2 = dAij(in, f, c);
			if (val2 > maxVal)
				maxVal = val2;
		}
		dAi(out,c) = maxVal;
	}
}

/* gets maximum value for each column*/
void ProgAngularAssignmentMag::meanByColumn(MultidimArray<double> &in,
		MultidimArray<double> &out) {

	out.resizeNoCopy(1, XSIZE(in));
	int f, c;
	double val, val2;
	int factor = YSIZE(in);
	for (c = 0; c < XSIZE(in); c++) {
		val = dAij(in, 0, c);
		for (f = 1; f < YSIZE(in); f++) {
			val2 = dAij(in, f, c);
			val += val2 / factor;
		}
		dAi(out,c) = val;
	}
}

/* gets maximum value for each row */
void ProgAngularAssignmentMag::maxByRow(MultidimArray<double> &in,
		MultidimArray<double> &out) {
	out.resizeNoCopy(1, YSIZE(in));
	int f, c;
	double maxVal, val2;
	for (f = 0; f < YSIZE(in); f++) {
		maxVal = dAij(in, f, 0);
		for (c = 1; c < XSIZE(in); c++) {
			val2 = dAij(in, f, c);
			if (val2 > maxVal)
				maxVal = val2;
		}
		dAi(out,f) = maxVal;
	}
}

/* gets maximum value for each row */
void ProgAngularAssignmentMag::meanByRow(MultidimArray<double> &in,
		MultidimArray<double> &out) {
	out.resizeNoCopy(1, YSIZE(in));
	int f, c;
	double val, val2;
	int factor = XSIZE(in);
	for (f = 0; f < YSIZE(in); f++) {
		val = dAij(in, f, 0);
		for (c = 1; c < XSIZE(in); c++) {
			val2 = dAij(in, f, c);
			val += val2 / factor;
		}
		dAi(out,f) = val;
	}
}

/*quadratic interpolation for location of peak in crossCorr vector*/
double quadInterp(/*const*/int idx, MultidimArray<double> &in) {
	int idxBef = idx - 1;
	int idxAft = idx + 1;
	if (idxBef < 0)
		std::cout << "Heeey!!\n";
	double oneBefore = dAi(in, idxBef);
	double oneAfter = dAi(in, idxAft);
	double centralOne = dAi(in, idx);
	double InterpIdx = idx
			- (( dAi(in,idx+1) - dAi(in, idx - 1))
					/ ( dAi(in,idx+1) + dAi(in, idx - 1) - 2 * dAi(in, idx)))
					/ 2.;
	return InterpIdx;
}

/* precompute Hann 2D window Ydim x Xdim*/
void ProgAngularAssignmentMag::computeHann() {
	float pi = 3.141592653;
	W.resizeNoCopy(Ydim, Xdim);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(W)
	{
		dAij(W,i,j) = 0.25 * (1 - cos(2 * pi * i / Xdim))
				* (1 - cos(2 * pi * j / Ydim));
	}
}

/*apply hann window to input image*/
void ProgAngularAssignmentMag::hannWindow(MultidimArray<double> &in) {
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(in)
	{
		dAij(in,i,j) *= dAij(W, i, j);
	}
}

/* precompute circular 2D window Ydim x Xdim*/
void ProgAngularAssignmentMag::computeCircular() {

	double Cf = (Ydim + (Ydim % 2)) / 2.0; // for odd/even images
	double Cc = (Xdim + (Xdim % 2)) / 2.0;
	int pixReduc = 1;
	double rad2 = (Cf - pixReduc) * (Cf - pixReduc);
	double val = 0;

	C.resizeNoCopy(Ydim, Xdim);
	C.initZeros(Ydim, Xdim);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(C)
	{
		val = (j - Cf) * (j - Cf) + (i - Cc) * (i - Cc);
		if (val < rad2)
			dAij(C,i,j) = 1.;
	}
}

/*apply circular window to input image*/
void ProgAngularAssignmentMag::circularWindow(MultidimArray<double> &in) {
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(in)
	{
		dAij(in,i,j) *= dAij(C, i, j);
	}
}

/* Only for 180 angles */
/* approach which selects only two locations of maximum peaks in ccvRot */
void ProgAngularAssignmentMag::rotCandidates3(MultidimArray<double> &in,
		std::vector<double> &cand, const size_t &size) {
	double max1 = -1000.;
	int idx1 = 0;
	double max2 = -1000.;
	int idx2 = 0;
	int i;
	int cont = 0;
	peaksFound = cont;

	for (i = 89; i < 272; ++i) { // only look for in range 90:-90
		// current value is a peak value?
		if ((dAi(in,size_t(i)) > dAi(in, size_t(i - 1)))
				&& (dAi(in,size_t(i)) > dAi(in, size_t(i + 1)))) {
			if ( dAi(in,i) > max1) {
				max2 = max1;
				idx2 = idx1;
				max1 = dAi(in, i);
				idx1 = i;
				cont += 1;
			} else if ( dAi(in,i) > max2 && dAi(in,i) != max1) {
				max2 = dAi(in, i);
				idx2 = i;
				cont += 1;
			}
		}
	}
	if (idx1 != 0) {
		int maxAccepted = 1;
		std::vector<int> temp;
		if (idx2 != 0) {
			maxAccepted = 2;
			temp.resize(maxAccepted);
			temp[0] = idx1;
			temp[1] = idx2;
		} else {
			temp.resize(maxAccepted);
			temp[0] = idx1;
		}

		//		// errase. Just checking what about with only one max
		//		int maxAccepted=1;
		//		std::vector<int> temp(maxAccepted);
		//		temp[0]=idx1;

		int tam = 2 * maxAccepted;
		peaksFound = tam;
		double interpIdx; // quadratic interpolated location of peak
		for (i = 0; i < maxAccepted; ++i) {
			interpIdx = quadInterp(temp[i], in);
			cand[i] = double(size) / 2. - interpIdx;
			cand[i + maxAccepted] =
					(cand[i] >= 0) ? cand[i] + 180. : cand[i] - 180.;
		}
	} else {
		peaksFound = 0;
	}
}

/* several candidates to best angle psi between images
 * then they are sorted and only keep maxAccepted
 * */
void ProgAngularAssignmentMag::rotCandidates(MultidimArray<double> &in,
		std::vector<double> &cand, const size_t &size) {
	int maxAccepted = 4;
	int maxNumOfPeaks = 90; // in 180 range is the max num of peaks, or even less
	std::vector<int> peakPos(maxNumOfPeaks, 0);
	std::vector<int> peakIdx(maxNumOfPeaks, 0);
	int cont = 0;
	peaksFound = cont;
	int i;
	for (i = 89/*1*/; i < 272/*size-1*/; ++i) { // check only the range -90:90
		if ((dAi(in,i) > dAi(in, i - 1)) && (dAi(in,i) > dAi(in, i + 1))) {
			peakIdx[cont] = cont; // for posterior ordering
			peakPos[cont] = i; // position of peak
			cont++;
			peaksFound = cont;
		}
	}

	maxAccepted = (peaksFound < maxAccepted) ? peaksFound : maxAccepted;

	if (cont) {
		// sort //todo check if its better to use partial_sort
		std::sort(peakIdx.begin(), peakIdx.end(),
				[&](int i, int j) {return dAi(in,i) > dAi(in,j);});
		//change for partial sort
		//		std::partial_sort(temp.begin(), temp.begin()+maxAccepted, temp.end(),
		//				[&](int i, int j){return dAi(in,i) > dAi(in,j); }); //

		int tam = 2 * maxAccepted; //
		peaksFound = tam;
		double interpIdx; // quadratic interpolated location of peak
		for (i = 0; i < maxAccepted; ++i) {
			interpIdx = quadInterp(peakPos[peakIdx[i]], in);
			cand[i] = double(size) / 2. - interpIdx;
			cand[i + maxAccepted] =
					(cand[i] >= 0) ? cand[i] + 180. : cand[i] - 180.;
		}
	} else {
		peaksFound = 0;
	}
}

/* instance of "delay axes" for assign rotation and translation candidates*/
void ProgAngularAssignmentMag::_delayAxes(const size_t &Ydim,
		const size_t &Xdim, const size_t &n_ang) {
	axRot.resize(1, 1, 1, n_ang);
	axTx.resize(1, 1, 1, Xdim);
	axTy.resize(1, 1, 1, Ydim);

	double M = double(n_ang - 1) / 2.;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axRot)
	{
		dAi(axRot,i) = ceil(M - i);
	}
	M = double(Xdim - 1) / 2.0;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axTx)
	{
		dAi(axTx,i) = ceil(M - i);
	}
	M = double(Ydim - 1) / 2.0;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(axTy)
	{
		dAi(axTy,i) = ceil(M - i);
	}
}

/* selection of best candidate to rotation and its corresponding shift
 * called at first loop in "coarse" searching
 * shifts are computed as maximum of CrossCorr vector
 * vector<double> cand contains candidates to relative rotation between images
 */
void ProgAngularAssignmentMag::bestCand(/*inputs*/
const MultidimArray<double> &MDaIn,
		const MultidimArray<std::complex<double> > &MDaInF,
		const MultidimArray<double> &MDaRef, std::vector<double> &cand,
		/*outputs*/
		double &bestPsi, double &shift_x, double &shift_y, double &bestCoeff) {
	bestPsi = 0;
	shift_x = 0.;
	shift_y = 0.;
	bestCoeff = 0.0;
	double rotVar = 0.0;
	double tempCoeff;
	double tx, ty;
	MultidimArray<double> MDaRefRot;
	MultidimArray<double> MDaRefShiftRot;
	MultidimArray<double> ccMatrixShift;
	MultidimArray<double> ccVectorTx;
	MultidimArray<double> ccVectorTy;
	MultidimArray<std::complex<double> > MDaRefRotF;

	MDaRefRot.setXmippOrigin();
	MDaRefShiftRot.setXmippOrigin();

	MultidimArray<double> MDaInShift;
	MultidimArray<double> MDaInShiftRot;
	MDaInShift.setXmippOrigin();
	MDaInShiftRot.setXmippOrigin();

	for (int i = 0; i < peaksFound; ++i) {
		rotVar = -1. * cand[i];  //negative, because is for reference rotation
		_applyRotation(MDaRef, rotVar, MDaRefRot); //rotation to reference image
		_applyFourierImage2(MDaRefRot, MDaRefRotF); //fourier
		ccMatrix(MDaInF, MDaRefRotF, ccMatrixShift); // cross-correlation matrix
		maxByColumn(ccMatrixShift, ccVectorTx); // ccvMatrix to ccVector
		getShift(ccVectorTx, tx, XSIZE(ccMatrixShift));
		tx = -1. * tx;
		maxByRow(ccMatrixShift, ccVectorTy); // ccvMatrix to ccVector
		getShift(ccVectorTy, ty, YSIZE(ccMatrixShift));
		ty = -1. * ty;

		if (std::abs(tx) > maxShift || std::abs(ty) > maxShift)
			continue;

		//apply transformation to experimental image
		double expTx, expTy, expPsi;
		expPsi = -rotVar;
		expTx = -tx;
		expTy = -ty;
		// applying in one transform
		_applyShiftAndRotation(MDaIn, expPsi, expTx, expTy, MDaInShiftRot);

		circularWindow(MDaInShiftRot); //circular masked MDaInRotShift
		// TODO compute another metric and check agreement between them in order to select better candidates
		// todo probably, candidate selection is good enough to at least have good candidates in top positions
		//      then in second loop could be useful to compute another metric, as a "second criteria" as in regularization
		//      could be useful to check information given by neighbors. e.g. low variance in alignment parameters in each neighborhood
		pearsonCorr(MDaRef, MDaInShiftRot, tempCoeff);  // Pearson
		//		normalized_cc(MDaRef, MDaInShiftRot, tempCoeff); // NCC
		//		imNormalized_cc(MDaRef, MDaInShiftRot, tempCoeff); // IMNCC
		//		imZNCC(MDaRef, MDaInShiftRot, tempCoeff); // IMZNCC
		// TODO COSS suggests use MDaRef as "mask"
		if (tempCoeff > bestCoeff) {
			bestCoeff = tempCoeff;
			shift_x = -expTx; //negative because in second loop,when used, this parameters are applied to mdaRef
			shift_y = -expTy;
			bestPsi = -expPsi;
		}
	}
}

/* same as bestCand but using 2 candidates for shift instead of 1 only
 */
void ProgAngularAssignmentMag::bestCand2(/*inputs*/
const MultidimArray<double> &MDaIn,
		const MultidimArray<std::complex<double> > &MDaInF,
		const MultidimArray<double> &MDaRef, std::vector<double> &cand,
		/*outputs*/
		double &bestPsi, double &shift_x, double &shift_y, double &bestCoeff) {
	bestPsi = 0;
	shift_x = 0.;
	shift_y = 0.;
	bestCoeff = 0.0;
	double rotVar = 0.0;
	double tempCoeff;
	double tx, ty;
	MultidimArray<double> MDaRefRot;
	MultidimArray<double> MDaRefShiftRot;
	MultidimArray<double> ccMatrixShift;
	MultidimArray<double> ccVectorTx;
	MultidimArray<double> ccVectorTy;
	MultidimArray<std::complex<double> > MDaRefRotF;

	MDaRefRot.setXmippOrigin();
	MDaRefShiftRot.setXmippOrigin();

	MultidimArray<double> MDaInShift;
	MultidimArray<double> MDaInShiftRot;
	MDaInShift.setXmippOrigin();
	MDaInShiftRot.setXmippOrigin();

	for (int i = 0; i < peaksFound; ++i) {
		rotVar = -1. * cand[i];  //negative, because is for reference rotation
		_applyRotation(MDaRef, rotVar, MDaRefRot); //rotation to reference image
		_applyFourierImage2(MDaRefRot, MDaRefRotF); //fourier
		ccMatrix(MDaInF, MDaRefRotF, ccMatrixShift); // cross-correlation matrix

		// testing 2 candidates for each psi
		std::vector<double> vTx(2, 0.);
		std::vector<double> vTy(2, 0.);
		maxByColumn(ccMatrixShift, ccVectorTx); // ccvMatrix to ccVector
		getShift2(ccVectorTx, vTx, XSIZE(ccMatrixShift));
		maxByRow(ccMatrixShift, ccVectorTy); // ccvMatrix to ccVector
		getShift2(ccVectorTy, vTy, YSIZE(ccMatrixShift));

		//      // todo In case this approach works properly, then set a condition to "continue"
		//		if (std::abs(vTx[0]) > maxShift || std::abs(vTy[0]) > maxShift ||
		//				std::abs(vTx[1]) > maxShift || std::abs(vTy[1]) > maxShift)
		//			continue;

		double expTx, expTy, expPsi;
		expPsi = -rotVar;
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				//apply transformation to experimental image
				expTx = vTx[j];
				expTy = vTy[k];
				// applying in one transform
				_applyShiftAndRotation(MDaIn, expPsi, expTx, expTy,
						MDaInShiftRot);
				circularWindow(MDaInShiftRot); //circular masked MDaInRotShift
				pearsonCorr(MDaRef, MDaInShiftRot, tempCoeff);  // Pearson
				if (tempCoeff > bestCoeff) {
					bestCoeff = tempCoeff;
					shift_x = -expTx; //negative because in second loop,when used, this parameters are applied to mdaRef
					shift_y = -expTy;
					bestPsi = -expPsi;
				}
			}
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
void ProgAngularAssignmentMag::newApplyGeometry(
		MultidimArray<double>& __restrict__ in,
		MultidimArray<double>& __restrict__ out, const double &a,
		const double &b, const double &c, const double &d, const double &tx,
		const double &ty) {

	int nFil = YSIZE(in);
	int nCol = XSIZE(in);

	double Cx = (nCol) / 2.0;
	double Cy = (nFil) / 2.0;

	// constants
	double k1 = b * Cy;
	double k2 = d * Cy;

	double e1 = Cx + tx;
	double e2 = Cy + ty;

	double d1 = e1 + k1; // Cx + tx + b * Cy
	double d2 = e2 + k2; // Cy + ty + d * Cy

	double g1 = e1 - k1; // Cx + tx - b * Cy
	double g2 = e2 - k2; // Cy + ty - d * Cy

	double x1, y1, p1, q1, p2, q2, p3, q3, p4, q4;

	int x, y, rx, ry;

	int lim_x1 = 0;
	int lim_x2 = nCol - lim_x1;
	int lim_y1 = 1;
	int lim_y2 = nFil - lim_y1;

	for (x = 0.; x < Cx; x++) {
		for (y = 0; y < Cy; y++) {
			x1 = a * double(x) + b * double(y);
			y1 = c * double(x) + d * double(y);

			// point 1 (x,y) // 4th
			p1 = x1 + e1;
			q1 = y1 + e2;
			rx = x + Cx;
			ry = y + Cy;
			if ((p1 > lim_x1) && (p1 < lim_x2) && (q1 > lim_y1)
					&& (q1 < lim_y2)) {
				dAij(out, ry, rx) = interpolate(in, q1, p1);
			}

			// point 2 (-x, -y + Cy) // 3th
			p2 = -x1 + d1;
			q2 = -y1 + d2;
			rx = -x + Cx;
			ry = -y + 2 * Cy;
			if ((p2 > lim_x1) && (p2 < lim_x2) && (q2 > lim_y1) && (q2 < lim_y2)
					&& (ry < lim_y2)) {
				dAij(out, ry, rx) = interpolate(in, q2, p2);
			}

			//point 3 (-x, -y) // 2nd
			p3 = -x1 + e1;
			q3 = -y1 + e2;
			rx = -x + Cx;
			ry = -y + Cy;
			if ((p3 > lim_x1) && (p3 < lim_x2) && (q3 > lim_y1)
					&& (q3 < lim_y2)) {
				dAij(out, ry, rx) = interpolate(in, q3, p3);
			}

			// point 4 (x, y-Cy) // 1st
			p4 = x1 + g1;
			q4 = y1 + g2;
			rx = x + Cx;
			ry = y;
			if ((p4 > lim_x1) && (p4 < lim_x2) && (q4 > lim_y1)
					&& (q4 < lim_y2)) {
				dAij(out, ry, rx) = interpolate(in, q4, p4);
			}
		}
	}

}

/* apply rotation */
void ProgAngularAssignmentMag::_applyRotation(
		const MultidimArray<double> &MDaRef, double &rot,
		MultidimArray<double> &MDaRefRot) {
	// Transform matrix
	Matrix2D<double> A(3, 3);
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

/* apply rotation */
void ProgAngularAssignmentMag::_applyRotation(MultidimArray<double> &MDaRef,
		double &rot, MultidimArray<double> &MDaRefRot) {
	// Transform matrix
	Matrix2D<double> A(3, 3);
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
void ProgAngularAssignmentMag::_applyShift(MultidimArray<double> &input,
		double &tx, double &ty, MultidimArray<double> &output) {
	// Transform matrix
	Matrix2D<double> A(3, 3);
	A.initIdentity();

	// Shift
	MAT_ELEM(A,0, 2) = tx;
	MAT_ELEM(A,1, 2) = ty;

	applyGeometry(LINEAR, output, input, A, IS_NOT_INV, DONT_WRAP);
}

/* apply traslation */
void ProgAngularAssignmentMag::_applyShift(const MultidimArray<double> &input,
		double &tx, double &ty, MultidimArray<double> &output) {
	// Transform matrix
	Matrix2D<double> A(3, 3);
	A.initIdentity();

	// Shift
	MAT_ELEM(A,0, 2) = tx;
	MAT_ELEM(A,1, 2) = ty;

	applyGeometry(LINEAR, output, input, A, IS_NOT_INV, DONT_WRAP);
}

/* finds shift as maximum of ccVector */
void ProgAngularAssignmentMag::getShift(MultidimArray<double> &ccVector,
		double &shift, const size_t &size) {
	double maxVal = -10.;
	int idx = 0;
	int i;
	int lb = int(size / 2 - maxShift);
	int hb = int(size / 2 + maxShift);
	for (i = 1; i < size - 1; ++i) { //i = lb; i < hb; ++i
		if (( dAi(ccVector,size_t(i)) > dAi(ccVector, size_t(i - 1)))
				&& ( dAi(ccVector,size_t(i)) > dAi(ccVector, size_t(i + 1)))) { // is this value a peak value?
			if (dAi(ccVector,i) > maxVal) { // is the biggest?
				maxVal = dAi(ccVector, i);
				idx = i;
			}
		}
	}

	if (idx) {
		// interpolate value
		double interpIdx;
		interpIdx = quadInterp(idx, ccVector);
		shift = double(size) / 2. - interpIdx;
	} else {
		shift = 0;
	}

}

/* finds shift as maximum of ccVector */
void ProgAngularAssignmentMag::getShift2(MultidimArray<double> &ccVector,
		std::vector<double> &cand, const size_t &size) {
	double max1 = -1000.;
	int idx1 = 0;
	double max2 = -1000.;
	int idx2 = 0;
	int i;

	for (i = 1; i < size - 1; ++i) { // only look for in range 90:-90
		// current value is a peak value?
		if ((dAi(ccVector,size_t(i)) > dAi(ccVector, size_t(i - 1)))
				&& (dAi(ccVector,size_t(i)) > dAi(ccVector, size_t(i + 1)))) {
			if ( dAi(ccVector,i) > max1) {
				max2 = max1;
				idx2 = idx1;
				max1 = dAi(ccVector, i);
				idx1 = i;
			} else if ( dAi(ccVector,i) > max2 && dAi(ccVector,i) != max1) {
				max2 = dAi(ccVector, i);
				idx2 = i;
			}
		}
	}
	if (idx1 != 0) {
		int maxAccepted = 1;
		std::vector<int> temp;
		if (idx2 != 0) {
			maxAccepted = 2;
			temp.resize(maxAccepted);
			temp[0] = idx1;
			temp[1] = idx2;
		} else {
			temp.resize(maxAccepted);
			temp[0] = idx1;
		}

		double interpIdx; // quadratic interpolated location of peak
		for (i = 0; i < maxAccepted; ++i) {
			interpIdx = quadInterp(temp[i], ccVector);
			cand[i] = double(size) / 2. - interpIdx;
		}
	} else {
		std::cout << "no peaks!\n";
	}
}

/* finds rot as maximum of ccVector for a region near the center */
void ProgAngularAssignmentMag::getRot(MultidimArray<double> &ccVector,
		double &rot, const size_t &size) {
	double maxVal = -10.;
	int idx;
	int i;
	//    int lb= int(size/2-5);
	//    int hb=int(size/2+5);
	int lb = 89;
	int hb = 270;
	for (i = lb; i < hb + 1; ++i) {
		if (dAi(ccVector,i) > maxVal) {
			maxVal = dAi(ccVector, i);
			idx = i;
		}
	}
	// interpolate value
	double interpIdx;
	interpIdx = quadInterp(idx, ccVector);
	//    rot = double( size - 1 )/2. - interpIdx;
	rot = double(size) / 2. - interpIdx;
}

/* Structural similarity SSIM index Coeff */
void ProgAngularAssignmentMag::ssimIndex(MultidimArray<double> &X,
		MultidimArray<double> &Y, double &coeff) {

	// covariance
	double X_m, Y_m, X_std, Y_std;
	double c1, c2, L;
	arithmetic_mean_and_stddev(X, X_m, X_std);
	arithmetic_mean_and_stddev(Y, Y_m, Y_std);

	double prod_mean = mean_of_products(X, Y);
	double covariace = prod_mean - (X_m * Y_m);

	L = 1;
	c1 = (0.01 * L) * (0.01 * L);
	c2 = (0.03 * L) * (0.03 * L); // estabilidad en división

	coeff = ((2 * X_m * Y_m + c1) * (2 * covariace + c2))
			/ ((X_m * X_m + Y_m * Y_m + c1)
					* (X_std * X_std + Y_std * Y_std + c2));
}

/* Structural similarity SSIM index Coeff */
void ProgAngularAssignmentMag::ssimIndex(const MultidimArray<double> &X,
		MultidimArray<double> &Y, double &coeff) {

	// covariance
	double X_m, Y_m, X_std, Y_std;
	double c1, c2, L;
	arithmetic_mean_and_stddev(X, X_m, X_std);
	arithmetic_mean_and_stddev(Y, Y_m, Y_std);

	double prod_mean = mean_of_products(X, Y);
	double covariace = prod_mean - (X_m * Y_m);

	L = 1;
	c1 = (0.01 * L) * (0.01 * L);
	c2 = (0.03 * L) * (0.03 * L); // stability in division

	coeff = ((2 * X_m * Y_m + c1) * (2 * covariace + c2))
			/ ((X_m * X_m + Y_m * Y_m + c1)
					* (X_std * X_std + Y_std * Y_std + c2));
}

/* apply rotation then shift */
void ProgAngularAssignmentMag::_applyRotationAndShift(
		const MultidimArray<double> &MDaRef, double &rot, double &tx,
		double &ty, MultidimArray<double> &MDaRefRot) {
	// Transform matrix
	Matrix2D<double> A(3, 3);
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

/* apply shift then rotation */
void ProgAngularAssignmentMag::_applyShiftAndRotation(
		const MultidimArray<double> &MDaRef, double &rot, double &tx,
		double &ty, MultidimArray<double> &MDaRefRot) {
	// Transform matrix
	Matrix2D<double> A(3, 3);
	A.initIdentity();
	double ang, cosine, sine;
	ang = DEG2RAD(rot);
	cosine = cos(ang);
	sine = sin(ang);

	// rotate in opposite direction
	double realTx = cosine * tx + sine * ty;
	double realTy = -sine * tx + cosine * ty;

	// rotation
	MAT_ELEM(A,0, 0) = cosine;
	MAT_ELEM(A,0, 1) = sine;
	MAT_ELEM(A,1, 0) = -sine;
	MAT_ELEM(A,1, 1) = cosine;

	// Shift
	MAT_ELEM(A,0, 2) = realTx;
	MAT_ELEM(A,1, 2) = realTy;

	applyGeometry(LINEAR, MDaRefRot, MDaRef, A, IS_NOT_INV, DONT_WRAP);

}

