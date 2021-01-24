



//  HERE IS THE LIKELIHOOD FOR THE COX REGRESSION
struct coxlikelihood{
	MatDoub simdata;
	coxlikelihood ();
	void Setcoxlikelihood( MatDoub &simdata );
	VecDoub operator()(const VecDoub_I x);

};

coxlikelihood::coxlikelihood() { }

void coxlikelihood::Setcoxlikelihood( MatDoub& simdata) {
	simdata = simdata;
}

VecDoub coxlikelihood::operator()(const VecDoub x) {
	int kkdn = simdata.nrows();
	VecDoub myreturn(4);
	int i;
	Doub thislambdai, thisbmkr, thistrt, thissurvtime, thissurvstat;
	for (i=0;i<4;i++) { myreturn[i]  = 0.0; }
	for (i=0;i<kkdn;i++) {
		thisbmkr = simdata[i][0];
		thistrt = simdata[i][1];
		thissurvtime = simdata[i][2];
		thissurvstat = simdata[i][3];
		thislambdai = exp(x[0] + x[1]*thisbmkr + x[2]*thistrt + x[3]*thisbmkr*thistrt);

		myreturn[0] = myreturn[0] - thissurvstat + thissurvtime * thislambdai;
		myreturn[1] = myreturn[1] - thissurvstat*thisbmkr + thissurvtime*thisbmkr*thislambdai;
		myreturn[2] = myreturn[2] - thissurvstat*thistrt + thissurvtime*thistrt*thislambdai;
		myreturn[3] = myreturn[3] - thissurvstat*thistrt*thisbmkr + thissurvtime*thistrt*thisbmkr*thislambdai;
	}
	return myreturn;
}
