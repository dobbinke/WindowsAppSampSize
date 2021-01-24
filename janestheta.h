
// HERE IS THE JANESTHETA STRUCT
struct JanesTheta {
public:
	MatDoub simdata;
	Doub t0;
	Doub beta0hatlogit, beta1hatlogit, beta2hatlogit, beta3hatlogit;
	Int kkdn;
	JanesTheta();
	void setJanesTheta ( MatDoub_I & , Doub & ,
	Doub & inputbeta0, Doub & inputbeta1, Doub & inputbeta2, Doub& inputbeta3);
	Doub Risk0 ( Doub thisbmkrval );
	Doub Risk1 ( Doub thisbmkrval );
	Doub Deltahat ( Doub thisbmkrval );
	Doub ProbDis1givenAis1andDeltaYover0 ( );
	Doub ProbDis1givenAis0andDeltaYover0 ( );
	Doub ProbDis1givenAis1andDeltaYunder0 ( );
	Doub ProbDis1givenAis0andDeltaYunder0 ( );
	Doub hatBsubneg( );
	Doub hatBsubpos( );
	Doub phatsubneg();
	Doub phatsubpos();
	Doub theta0hat ();
	Doub theta1hat ();
	MatDoub PercentileBootstrap (VecDoub_I percentiles, int nboots);
	MatDoub DoTheRefinedBootstrap (int nboots);
};

JanesTheta::JanesTheta() {}

void JanesTheta::setJanesTheta(MatDoub_I &inputdataset, Doub & inputt0,
	Doub & inputbeta0, Doub & inputbeta1, Doub & inputbeta2, Doub& inputbeta3)
{
    MatDoub simdata = inputdataset;
    kkdn = inputdataset.nrows();
	t0 = inputt0;
	beta0hatlogit = inputbeta0;
	beta1hatlogit = inputbeta1;
	beta2hatlogit = inputbeta2;
	beta3hatlogit = inputbeta3;
}

Doub JanesTheta::Risk0( Doub thisbmkrval ) {
	Doub thisreturn;
//	thisreturn = 1.0 - exp( -exp( beta0hatlogit + beta1hatlogit * thisbmkrval)); //Cox regression
	Doub myl1 = beta0hatlogit + beta1hatlogit * thisbmkrval;
	thisreturn = exp( myl1 )/(1+exp( myl1 )); // regression

	return thisreturn;
}

Doub JanesTheta::Risk1( Doub thisbmkrval ) {
	Doub thisreturn;

//	thisreturn = 1.0 - exp( -exp( beta0hatlogit + beta1hatlogit * thisbmkrval + beta2hatlogit + beta3hatlogit * thisbmkrval));
	Doub myl1 = beta0hatlogit + beta1hatlogit * thisbmkrval + beta2hatlogit + beta3hatlogit * thisbmkrval;
	thisreturn = exp( myl1 )/(1+exp( myl1 )); //logistic regression

	return thisreturn;
}

Doub JanesTheta::Deltahat( Doub thisbmkrval ) {
	Doub thisreturn;
	thisreturn = Risk0(thisbmkrval) - Risk1(thisbmkrval);
	return thisreturn;
}

Doub JanesTheta::ProbDis1givenAis1andDeltaYover0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;


	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat > 0);
		Bool condition3 = (thistrt == 1);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = mynumerator/mydenominator;

	return myreturn;
}

Doub JanesTheta::ProbDis1givenAis0andDeltaYover0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;


	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat > 0);
		Bool condition3 = (thistrt == 0);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = mynumerator/mydenominator;

	return myreturn;
}


Doub JanesTheta::hatBsubneg() {
	Doub myreturn;
	myreturn =  ProbDis1givenAis1andDeltaYunder0() - ProbDis1givenAis0andDeltaYunder0();
	return myreturn;
}

Doub JanesTheta::hatBsubpos() {
	Doub myreturn;
	myreturn =  ProbDis1givenAis0andDeltaYover0() - ProbDis1givenAis1andDeltaYover0();
	return myreturn;
}


Doub JanesTheta::phatsubneg( ) {
	Doub myreturn;
	Doub mynum= 0, mydenom = 0;
	Doub thistime, thissurvstat, thisbmkrval, thiscensored, thisnotcensored;

	for (int i=0;i<kkdn;i++) {
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisbmkrval = simdata[i][0];

		thiscensored = (thistime<5)*(thissurvstat==0);
		thisnotcensored = 1.-thiscensored;

		mydenom = thisnotcensored + mydenom;
		mynum = mynum + thisnotcensored*(Deltahat(thisbmkrval)<0);
	}
	myreturn = mynum/mydenom;

	return myreturn;
}

Doub JanesTheta::phatsubpos( ) {
	Doub myreturn;
	myreturn = 1. - phatsubneg();
	return myreturn;
}

Doub JanesTheta::theta1hat( ) {
	Doub myreturn;
	myreturn = hatBsubneg() * phatsubneg();
	// myreturn = max(myreturn,0.0);
	return myreturn;
}

Doub JanesTheta::theta0hat( ) {
	Doub myreturn;
	myreturn = hatBsubpos() * phatsubpos();
	// myreturn = max(myreturn,0.0);
	return myreturn;
}


Doub JanesTheta::ProbDis1givenAis1andDeltaYunder0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;
	kkdn = simdata.nrows();

	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat < 0);
		Bool condition3 = (thistrt == 1);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = mynumerator/mydenominator;

	return myreturn;
}

Doub JanesTheta::ProbDis1givenAis0andDeltaYunder0( ) {
	Doub myreturn;

	double mynumerator = 0.0, mydenominator = 0.0;
	Doub thisdeltahat, thisbmkrval, thistime, thistrt, thissurvstat, thisnotcensored;
	Doub thisnumsubset, thisdenomsubset;
	kkdn = simdata.nrows();

	for (int i=0;i<kkdn;i++) {
		thisbmkrval = simdata[i][0];
		thistrt = simdata[i][1];
		thistime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thisdeltahat = Deltahat(thisbmkrval);

		Bool condition1 = ((thistime > t0) || (thissurvstat == 1));
		Bool condition2 = (thisdeltahat < 0);
		Bool condition3 = (thistrt == 0);
		if (condition1 && condition2 && condition3) {
			mydenominator = mydenominator + 1.0;
		}
		Bool condition4 = (thistime < t0);
		if (condition1 && condition2 && condition3 && condition4) {
			mynumerator = mynumerator + 1.0;
		}
	}
	myreturn = double(mynumerator) / double(max(mydenominator,0.01));

	return myreturn;
}



MatDoub JanesTheta::DoTheRefinedBootstrap (int nboots) {
	JanesTheta resubjt, bootjt, bootjtresamp;
	Doub A0, A1, B0, B1, C0, C1, D0, D1;
	MatDoub myreturn(2,nboots);
	// first row is theta0s, second row is theta1s
	coxlikelihood tempcoxlikelihood;
	MatDoub BootstrapSimData(kkdn,4);
	int i, thisboot, thisindex;
	VecDoub betaMLE(4);
	Bool tempbool;

	// RESUBSTITUTION ESTIMATES
	tempcoxlikelihood.simdata = simdata;
	betaMLE[0] = 0; betaMLE[1] = 0; betaMLE[2] = 0; betaMLE[3] = 0;
	newt(betaMLE,tempbool,tempcoxlikelihood);
	resubjt.simdata = simdata; // the full dataset is used
	resubjt.setJanesTheta(simdata,t0,betaMLE[0],betaMLE[1],betaMLE[2],betaMLE[3]);
	   // the beta estimates are taken from the bootstrap fit
	D0 = resubjt.theta0hat();
	D1 = resubjt.theta1hat();

    // END OF RESUBSTITUTEION ESTIMATES


	for (thisboot=0;thisboot<nboots;thisboot++) {

// 1) Draw a bootstrap sample and develop a model.
		// Create bootstrap dataset
		for (i=0;i<kkdn;i++) {
			thisindex = rand() % kkdn;
			for (int j=0;j<4;j++) {
				BootstrapSimData[i][j] = simdata[thisindex][j];
			}
		}

		// Maximum likelihood Estimation
		tempcoxlikelihood.simdata = BootstrapSimData;
		betaMLE[0] = 0; betaMLE[1] = 0; betaMLE[2] = 0; betaMLE[3] = 0;
		newt(betaMLE,tempbool,tempcoxlikelihood);

// 2) Apply the model to all n samples in the dataset to estimate theta0 and theta1,
		// call these A0 and A1.
		bootjt.simdata = simdata; // the full dataset is used
		bootjt.setJanesTheta(simdata,t0,betaMLE[0],betaMLE[1],betaMLE[2],betaMLE[3]);
		   // the beta estimates are taken from the bootstrap fit
		A0 = bootjt.theta0hat();
		A1 = bootjt.theta1hat();

// 3) Apply the model just to the samples used to train the model (resubstitution error),
		//  call it B0 and B1.
		bootjtresamp.simdata = BootstrapSimData;
		bootjtresamp.setJanesTheta(BootstrapSimData,t0,betaMLE[0],betaMLE[1],betaMLE[2],betaMLE[3]);
		B0 = bootjtresamp.theta0hat();
		B1 = bootjtresamp.theta1hat();

// 4) The overoptimism bias is C = A-B.
		C0 = A0 - B0;
		C1 = A1 - B1;

// 5) Now develop a model on the full dataset and estimate the error with resubstitution, D0 and D1.

		// THIS IS ABOVE

// 6) The final refined bootstrap estimate is D+C.

		myreturn[0][thisboot] = D0+C0;
		myreturn[1][thisboot] = D1+C1;
	}  // END OF BOOTSTRAP LOOP

	return myreturn;
}


MatDoub JanesTheta::PercentileBootstrap ( VecDoub_I percentiles, int nboots ) {
	MatDoub myreturn(2,2);
	MatDoub BSresults(2,nboots);
	BSresults = DoTheRefinedBootstrap(nboots);
	VecDoub tempvec0(nboots);
	VecDoub tempvec1(nboots);
	int i;

	double thisfullsampleTheta0 = theta0hat();
    double thisfullsampleTheta1 = theta1hat();

	for (i=0;i<nboots;i++) {
		tempvec0[i] = BSresults[0][i];
		tempvec1[i] = BSresults[1][i];
	}

    // percentile bootstrap code;
	int rank1 = round(0.025 * nboots);
	int rank2 = round(0.975 * nboots);

	myreturn[0][0] = select(rank1, tempvec0);
	myreturn[0][1] = select(rank2, tempvec0);

	myreturn[1][0] = select(rank1, tempvec1);
	myreturn[1][1] = select(rank2, tempvec1);

	return  myreturn;
}
