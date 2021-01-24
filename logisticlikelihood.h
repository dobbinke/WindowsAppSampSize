

struct logisticlikelihood{
	MatDoub simdata;
	Doub t0;
	logisticlikelihood ();
	void Setlogisticlikelihood( MatDoub &simdata, Doub &t0 );
	VecDoub operator()(const VecDoub_I x);
};

logisticlikelihood::logisticlikelihood() { }

void logisticlikelihood::Setlogisticlikelihood( MatDoub& simdata, Doub &t0) {
	simdata = simdata;
	t0 = t0;
}

VecDoub logisticlikelihood::operator()(const VecDoub x) {
	int kkdn = simdata.nrows();
	VecDoub myreturn(4);
	int i;
	Doub thispii, thisli, thisbmkr, thistrt, thissurvtime, thissurvstat, thisBinaryResp;
	for (i=0;i<4;i++) { myreturn[i]  = 0.0; }
	for (i=0;i<kkdn;i++) {
		thisbmkr = simdata[i][0];
		thistrt = simdata[i][1];
		thissurvtime = simdata[i][2];
		thissurvstat = simdata[i][3];

		thisBinaryResp = 99;  // indicator of missing data
		if (thissurvtime < t0 && thissurvstat == 1) { thisBinaryResp = 1; } // bad outcome group
		if (thissurvtime >= t0) { thisBinaryResp = 0; } // good outcome group

		thisli = x[0] + x[1]*thisbmkr + x[2]*thistrt + x[3]*thisbmkr*thistrt;
		thispii = exp(thisli)/(1+exp(thisli));

		if (thisBinaryResp != 99) {
			myreturn[0] = myreturn[0] - thisBinaryResp + thispii;
			myreturn[1] = myreturn[1] - thisBinaryResp * thisbmkr + thisbmkr*thispii;
			myreturn[2] = myreturn[2] - thisBinaryResp * thistrt + thistrt*thispii;
			myreturn[3] = myreturn[3] - thisBinaryResp * thistrt * thisbmkr + thistrt*thisbmkr*thispii;
		}
	}
	return myreturn;
}
