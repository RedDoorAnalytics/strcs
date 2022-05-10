************************************************************************************	
* MATA CODE

version 1.82

local SS 	string scalar
local RS	real scalar
local NM	numeric matrix
local TM	transmorphic matrix
mata:

struct strcs_main {
	real colvector t, lnt, t0, d, xb, rcs, lowt, hight, lowerb, upperb, includefirstint, includesecondint, includethirdint, tvcoffsetsum, bhazard, bob
	`RS' delentry, Nobs, orthog, hastvc, hascovs, Nnodes, Ntvc, Nrcsterms, likeconstant, bhtime, rsmodel, nocons, Ncovs
	`SS' tvclist, tvcofflist, covslist
	`NM' tvcvars, weights
	`TM' invRmatrix, betapos, nodes, rcsfirstknot, rcslastknot, covsmat, rcst0knot
	real matrix cumhazardmid, basercs, basedrcs, cumhazard
}

void strcs_setup(`SS' temp)
{
	struct strcs_main scalar PS
	pointer scalar p
	rmexternal(temp)	
	p = crexternal(temp)

	if(st_local("verbose") != "") printf("{text:Setting up}\n")
	touse = st_local("touse")
	PS.hastvc = st_local("tvc") != ""
	PS.delentry = strtoreal(st_local("del_entry"))
	PS.Nobs = strtoreal(st_local("N"))
	PS.hascovs = st_local("varlist") != ""
	PS.orthog = st_local("orthog") != ""
	PS.bhtime = st_local("bhtime") != ""
	PS.rsmodel = st_local("bhazard") != ""
	PS.nocons = st_local("constant") != ""
	PS.basercs=st_data(.,tokens("__s*"),touse)
	PS.basedrcs=st_data(.,tokens("__d_s*"),touse)
	if(!PS.nocons) {
		PS.basercs= PS.basercs,J(PS.Nobs,1,0)
		PS.basedrcs=PS.basedrcs,J(PS.Nobs,1,0)
	}

	
	df = strtoreal(st_local("df"))	
	
	if (!PS.rsmodel) PS.likeconstant = strtoreal(st_local("likeconstant"))
	if (PS.hastvc) {
		PS.tvclist = tokens(st_local("tvc"))
		PS.Ntvc = cols(PS.tvclist)
	}
// Rmatrix
	if (PS.orthog) {
		baseterms = df + 1 - PS.nocons
		PS.invRmatrix = asarray_create()
		asarray(PS.invRmatrix,"_baseline",qrinv(st_matrix(st_local("R_bh"))[1..baseterms,1..baseterms])')
		if (PS.hastvc) {
			for(i=1;i<=PS.Ntvc;i++) {
				asarray(PS.invRmatrix,PS.tvclist[i],qrinv(st_matrix(st_local("R_"+PS.tvclist[i]))'))
			}
		}
	}

// Data
	PS.t = st_data(.,"_t",touse)
	PS.d = st_data(.,"_d",touse)
	PS.lowt = st_data(.,st_local("lowt"),touse)
	PS.hight = st_data(.,st_local("hight"),touse)
	PS.lowerb = st_data(.,st_local("lowerb"),touse)
	PS.upperb = st_data(.,st_local("upperb"),touse)
	if (!PS.rsmodel) PS.lnt = st_data(.,st_local("lnt"),touse)
	if (PS.rsmodel) PS.bhazard = st_data(.,st_local("bhazard"),touse)
	if (PS.delentry) PS.t0 = st_data(.,"_t0",touse)
	if (PS.hascovs) {
		PS.covslist = tokens(st_local("varlist")) 
		PS.Ncovs = cols(PS.covslist)
		PS.covsmat=st_data(.,PS.covslist,touse)
	}
	


	if (PS.hastvc) {
		PS.tvcvars = st_data(.,PS.tvclist,touse)
		if (st_local("tvcoffset") != "") PS.tvcoffsetsum= st_data(.,st_local("tvcoffsetsum"),touse) 
		else PS.tvcoffsetsum = 0
	}
	PS.includefirstint = st_data(.,st_local("includefirstint"),touse)
	PS.includesecondint = st_data(.,st_local("includesecondint"),touse)
	PS.includethirdint = st_data(.,st_local("includethirdint"),touse)
	

// nodes and weights
	PS.weights = J(PS.Nobs,1,st_matrix(st_local("kweights"))):*(PS.upperb :- PS.lowerb):/2
	PS.Nnodes = strtoreal(st_local("nodes"))
	knodes = st_matrix(st_local("knodes"))
	if (PS.hastvc) {
		tvc_df = J(1,PS.Ntvc,.)
		rcstvcname = J(1,PS.Ntvc,"") 
		drcstvcname = J(1,PS.Ntvc,"") 
		for(j=1;j<=PS.Ntvc;j++) {
			tvc_df[1,j] = strtoreal(st_local("tvc_"+PS.tvclist[j]+"_df"))
			rcstvcname[1,j] = "__rcs"+PS.tvclist[j] 
			drcstvcname[1,j] = "__drcs"+PS.tvclist[j] 
		}
	}
	
	PS.cumhazardmid = J(PS.Nobs,PS.Nnodes,0)	
	PS.cumhazard = J(PS.Nobs,PS.Nnodes,0)
	
	if(PS.nocons) cons=J(PS.Nobs,0,.)
	else cons = J(PS.Nobs,1,1)
	PS.nodes = asarray_create("real",1)
	PS.bhtime = st_local("bhtime") != ""
	
	st_view(tmpnode ,.,"__tmpnode",touse)
	st_view(tmpnodetvc,.,"__tmpnodetvc",touse)
	st_view(tmptvcgen,.,"__tmptvcgen",touse)
	stata("capture drop __rcs*")
	
	for(i=1;i<=PS.Nnodes;i++) {	
		tmpnode_i = (0.5:*(PS.upperb :- PS.lowerb):*knodes[1,i] :+ 0.5:*(PS.upperb :+ PS.lowerb)) 

		

		
		if (PS.bhtime) tmpnode[,] = tmpnode_i
		else tmpnode[,] = ln(tmpnode_i)
		stata("qui rcsgen __tmpnode if " + touse +", gen(__rcs) knots("+st_local("bhknots")+") "+st_local("R_bh_opt")+ " " +st_local("reverse")) 
		if (PS.hastvc) {
			for(j=1;j<=PS.Ntvc;j++) {
				
				if(st_local("tvc_"+PS.tvclist[j]+"_offset") != "") {
					st_view(offset,.,st_local("tvc_"+PS.tvclist[j]+"_offset"),touse)
				}
				else offset = 0
				
				if(PS.bhtime) tmpnodetvc[,] = tmpnode_i :+ offset
				else tmpnodetvc[,] = ln(tmpnode_i) :+ offset

				stata("qui rcsgen __tmpnodetvc if " + touse +", gen(__rcs"+PS.tvclist[j]+") knots("+st_local("tvcknots_"+PS.tvclist[j])+") " + st_local("R_tvc_"+PS.tvclist[j]+"_opt")+ " " +st_local("reverse")) 
				for(k=1;k<=tvc_df[1,j];k++) {
					stata("qui replace "+rcstvcname[1,j]+strofreal(k)+" =" +rcstvcname[j]+strofreal(k)+"*"+PS.tvclist[j]+" if "+touse)
				}
			}
		}
		asarray(PS.nodes,i,(st_data(.,tokens("__rcs*"),touse),cons))
		stata("capture drop __rcs* ")
	}
	stata("capture drop __tmpnode*")

//	rcs values at t0	
if(PS.delentry) {	
	PS.rcst0knot = asarray_create("real",1) // 1 - splines, 2 - derivatives
	
	if (PS.bhtime) stata("qui rcsgen " + st_local("t0") + " if " + touse +", gen(__rcs) dgen(__drcs) knots("+st_local("bhknots")+") "+st_local("R_bh_opt")+ " " +st_local("reverse"))
	else stata("qui rcsgen " + st_local("ln_t0") + " if " + touse +", gen(__rcs) dgen(__drcs) knots("+st_local("bhknots")+") "+st_local("R_bh_opt")+ " " +st_local("reverse"))
	if (PS.hastvc) {
		st_view(tmptvcgen,.,"__tmptvcgen",touse)
		for(j=1;j<=PS.Ntvc;j++) {
			if(st_local("tvc_"+PS.tvclist[j]+"_offset") != "") st_view(offset,.,st_local("tvc_"+PS.tvclist[j]+"_offset"),touse)
			else offset = 0	

			if (PS.bhtime) tmptvcgen[,] = PS.t0 :+ offset
			else tmptvcgen[,] = ln(PS.t0) :+ offset

			stata("qui rcsgen __tmptvcgen if " + touse +", gen(__rcs"+PS.tvclist[j]+") dgen(__drcs"+PS.tvclist[j]+") knots("+st_local("tvcknots_"+PS.tvclist[j])+") " + st_local("R_tvc_"+PS.tvclist[j]+"_opt")+ " " +st_local("reverse")) 
			for(k=1;k<=tvc_df[1,j];k++) {
				stata("qui replace "+rcstvcname[1,j]+strofreal(k)+" =" +rcstvcname[j]+strofreal(k)+"*"+PS.tvclist[j]+" if "+touse)
				stata("qui replace "+drcstvcname[1,j]+strofreal(k)+" =" +drcstvcname[j]+strofreal(k)+"*"+PS.tvclist[j]+" if "+touse)
		}
		}
	}
	asarray(PS.rcst0knot,1,(st_data(.,"__rcs*",touse),cons))

	if(!PS.nocons) {
		asarray(PS.rcst0knot,2,(st_data(.,"__drcs*",touse),J(PS.Nobs,1,0)))
	}
	else {
		asarray(PS.rcst0knot,2,(st_data(.,"__drcs*",touse)))
	}
	
	stata("capture drop __rcs*")
	stata("capture drop __drcs*")
			
	}

	
//	rcs values at first knot 
	PS.rcsfirstknot = asarray_create("real",1) // 1 - splines, 2 - derivatives
	
	
	if (PS.bhtime) 	stata("qui rcsgen " + st_local("lowt") + " if " + touse +", gen(__rcs) dgen(__drcs) knots("+st_local("bhknots")+") "+st_local("R_bh_opt")+ " " +st_local("reverse"))
	else stata("qui rcsgen " + st_local("ln_lowt") + " if " + touse +", gen(__rcs) dgen(__drcs) knots("+st_local("bhknots")+") "+st_local("R_bh_opt")+ " " +st_local("reverse"))
	
	if (PS.hastvc) {
		st_view(tmptvcgen,.,"__tmptvcgen",touse)
		for(j=1;j<=PS.Ntvc;j++) {
			if(st_local("tvc_"+PS.tvclist[j]+"_offset") != "") st_view(offset,.,st_local("tvc_"+PS.tvclist[j]+"_offset"),touse)
			else offset = 0	

			if (PS.bhtime) tmptvcgen[,] = PS.lowt :+ offset
			else tmptvcgen[,] = ln(PS.lowt) :+ offset

			stata("qui rcsgen __tmptvcgen if " + touse +", gen(__rcs"+PS.tvclist[j]+") dgen(__drcs"+PS.tvclist[j]+") knots("+st_local("tvcknots_"+PS.tvclist[j])+") " + st_local("R_tvc_"+PS.tvclist[j]+"_opt")+ " " +st_local("reverse")) 
			for(k=1;k<=tvc_df[1,j];k++) {
				stata("qui replace "+rcstvcname[1,j]+strofreal(k)+" =" +rcstvcname[j]+strofreal(k)+"*"+PS.tvclist[j]+" if "+touse)
				stata("qui replace "+drcstvcname[1,j]+strofreal(k)+" =" +drcstvcname[j]+strofreal(k)+"*"+PS.tvclist[j]+" if "+touse)
		}
		}
	}
	PD = asarray(PS.rcsfirstknot,1,(st_data(.,"__rcs*",touse),cons))

	if(!PS.nocons) {
		asarray(PS.rcsfirstknot,2,(st_data(.,"__drcs*",touse),J(PS.Nobs,1,0)))
	}
	else {
		asarray(PS.rcsfirstknot,2,(st_data(.,"__drcs*",touse)))
	}
	
	stata("capture drop __rcs*")
	stata("capture drop __drcs*")
	
//	rcs values at last knot 
	PS.rcslastknot = asarray_create("real",1) // 1 - splines, 2 - derivatives
	
	if (PS.bhtime) 	stata("qui rcsgen " + st_local("hight") + " if " + touse +", gen(__rcs) dgen(__drcs) knots("+st_local("bhknots")+") "+st_local("R_bh_opt")+ " " +st_local("reverse"))
	else stata("qui rcsgen " + st_local("ln_hight") + " if " + touse +", gen(__rcs) dgen(__drcs) knots("+st_local("bhknots")+") "+st_local("R_bh_opt")+ " " +st_local("reverse"))

	if (PS.hastvc) {
		st_view(tmptvcgen,.,"__tmptvcgen",touse)
		for(j=1;j<=PS.Ntvc;j++) {
			if(st_local("tvc_"+PS.tvclist[j]+"_offset") != "") st_view(offset,.,st_local("tvc_"+PS.tvclist[j]+"_offset"),touse)
			else offset = 0

			if (PS.bhtime) tmptvcgen[,] = PS.hight :+ offset
			else tmptvcgen[,] = ln(PS.hight) :+ offset

			stata("qui rcsgen __tmptvcgen if " + touse +", gen(__rcs"+PS.tvclist[j]+") dgen(__drcs"+PS.tvclist[j]+") knots("+st_local("tvcknots_"+PS.tvclist[j])+") " + st_local("R_tvc_"+PS.tvclist[j]+"_opt")+ " " +st_local("reverse")) 
			for(k=1;k<=tvc_df[1,j];k++) {
				stata("qui replace "+rcstvcname[1,j]+strofreal(k)+" =" +rcstvcname[j]+strofreal(k)+"*"+PS.tvclist[j]+" if "+touse)
				stata("qui replace "+drcstvcname[1,j]+strofreal(k)+" =" +drcstvcname[j]+strofreal(k)+"*"+PS.tvclist[j]+" if "+touse)
		}
		}
	}
	asarray(PS.rcslastknot,1,(st_data(.,"__rcs*",touse),cons))
	
	if(!PS.nocons) {
			asarray(PS.rcslastknot,2,(st_data(.,"__drcs*",touse),J(PS.Nobs,1,0)))
	}
	else {
		asarray(PS.rcslastknot,2,(st_data(.,"__drcs*",touse)))

	}
	
	stata("capture drop __rcs*")
	stata("capture drop __drcs*")
	stata("capture drop __tmptvcgen")
	
	// equations
	PS.Nrcsterms = strtoreal(st_local("Nrcsterms"))
	PS.betapos = asarray_create()

	if(!PS.nocons) addcons = PS.Nrcsterms
	else addcons = J(1,0,.)


	
	asarray(PS.betapos,"_baseline",(range(1,df,1)',addcons))	
	if (PS.hastvc) {
		tmppos  = df + 1
		for(i=1;i<=PS.Ntvc;i++) {
			dftvc = strtoreal(st_local("tvc_"+PS.tvclist[i]+"_df"))
			asarray(PS.betapos,PS.tvclist[i],(range(tmppos,(tmppos + dftvc - 1),1)',PS.Nrcsterms))	
		}
	}

	if (PS.hascovs) PS.xb = J(PS.Nobs,1,.)
	PS.rcs = J(PS.Nobs,1,.)
	
	if(st_local("verbose") != "") display("{text:Fitting Model}")
// Done	
swap((*p), PS)
}



function strcs_gf(transmorphic M, todo, b, lnf, S, H)
{
	pointer(struct strcs_main scalar) scalar PS
	PS = &moptimize_util_userinfo(M,1)
	
	if ((*PS).hascovs) {
		PS->xb = moptimize_util_xb(M,b,1)
		PS->rcs = moptimize_util_xb(M,b,2)
		rcsbeta = b[|moptimize_util_eq_indices(M, 2)|] 
		xb = b[|moptimize_util_eq_indices(M, 1)|] 
	
	}
	else {
		PS->rcs = moptimize_util_xb(M,b,1)
		PS->xb = 0
		rcsbeta = b[|moptimize_util_eq_indices(M, 1)|] 
		xb = 0

	}


// intercept and gradient before first knot
	b0 =  (asarray((*PS).rcsfirstknot,1) * rcsbeta') :+ (*PS).xb
	b1 = (asarray((*PS).rcsfirstknot,2) * rcsbeta') 
	if ((*PS).bhtime) b0 = b0 :- b1:*(*PS).lowt
	else b0 = b0 :- b1:*ln((*PS).lowt)	


// intercept and gradient after last knot
	b0_last =  (asarray((*PS).rcslastknot,1) * rcsbeta') :+ (*PS).xb
	b1_last = (asarray((*PS).rcslastknot,2) * rcsbeta') 
	
	if ((*PS).bhtime) b0_last = b0_last :- b1_last:*(*PS).hight
	else b0_last = b0_last :- b1_last:*ln((*PS).hight)	
	

// Up to first knot (analytic)- 
	if ((*PS).bhtime) {
		if (!(*PS).delentry) cumhazard = (exp(b0):/b1:*(exp(b1:*(*PS).lowt) :- 1)):*((*PS).includefirstint)
		else cumhazard = (exp(b0):/b1:*(exp(b1:*(*PS).lowt) :- exp(b1:*(*PS).t0))):*((*PS).includefirstint)
	}	
	else {
		if (!(*PS).delentry) cumhazard = (exp(b0):/ (b1:+1):*((*PS).lowt:^(b1:+1))):*((*PS).includefirstint)
		else cumhazard = (exp(b0):/ (b1:+1):*((*PS).lowt:^(b1:+1) :- (*PS).t0:^(b1:+1))):*((*PS).includefirstint)
	}

	
// Middle part (numerical)	

	for(i=1;i<=(*PS).Nnodes;i++) {
		PS-> cumhazardmid[,i] = ((*PS).weights[,i]:*exp((asarray((*PS).nodes, i)*rcsbeta' :+ (*PS).xb))) :*(*PS).includesecondint 
			PS-> cumhazard[,i] = ((*PS).weights[,i]:*exp((asarray((*PS).nodes, i)*rcsbeta' :+ (*PS).xb))) //:*(*PS).includesecondint
	}
	

	cumhazardmid = quadrowsum((*PS).cumhazardmid,1)
	cumhazard = cumhazard :+ cumhazardmid


// After last knot (analytical)
	if ((*PS).bhtime) cumhazard = cumhazard :+  (exp(b0_last):/b1_last:*(exp(b1_last:*(*PS).t) :- exp(b1_last:*(*PS).hight))):*((*PS).includethirdint)
	else cumhazard = cumhazard :+(exp(b0_last):/ (b1_last:+1):*((*PS).t:^(b1_last:+1) :- ((*PS).hight):^(b1_last:+1))):*((*PS).includethirdint)
		
// likelihood	
	if ((*PS).rsmodel) lnf = (*PS).d :*ln((*PS).bhazard :+ exp((*PS).rcs :+ (*PS).xb)) :- cumhazard
	else lnf = (*PS).d :*((*PS).rcs :+ (*PS).xb :+ (*PS).lnt) :- cumhazard  
	if (todo==0) return

	
	
	
////////////////////////////////////////////////////////////////////////////////////////////////	
//SCORE EQUATION
// create functions for the analytical part of the score equation before first knot
	u = (exp(b0))	
	if(!(*PS).bhtime) v = (1:/(b1 :+ 1))
	else v = (1:/b1)
	w=(b1:+1)
	
// create functions for analytical part of the score equation after the last knot	
	u_last = (exp(b0_last))
	if(!(*PS).bhtime) v_last = (1:/(b1_last :+ 1))
	else v_last = (1:/b1_last)
	w_last=(b1_last:+1)



// Up to first knot (analytic)
	if ((*PS).bhtime) {
		if (!(*PS).delentry) {
			cumhazardscoremin = (u:*v:*exp(b1:*(*PS).lowt):*asarray((*PS).rcsfirstknot,2):*(*PS).lowt ///
				:- u:*(v:^2):*asarray((*PS).rcsfirstknot,2):*exp(b1:*(*PS).lowt) ///
				:+ u:*(asarray((*PS).rcsfirstknot,1):-asarray((*PS).rcsfirstknot,2):*(*PS).lowt):*v:*exp(b1:*(*PS).lowt) ///
				:- u:*(asarray((*PS).rcsfirstknot,1) :- asarray((*PS).rcsfirstknot,2):*(*PS).lowt):*v ///
				:+ u:*(v:^2):*asarray((*PS).rcsfirstknot,2)):*((*PS).includefirstint)
		}
	else {
		cumhazardscoremin = ((u:*v:*exp(b1:*(*PS).lowt):*asarray((*PS).rcsfirstknot,2):*(*PS).lowt ///
			:-u:*(v:^2):*asarray((*PS).rcsfirstknot,2):*exp(b1:*(*PS).lowt) ///
			:+ u:*(asarray((*PS).rcsfirstknot,1) :- asarray((*PS).rcsfirstknot,2):*(*PS).lowt):*v:*exp(b1:*((*PS).lowt))) ///
			:- (u:*v:*exp(b1:*(*PS).t0):*asarray((*PS).rcst0knot,2):*(*PS).t0 ///
			:-u:*(v:^2):*asarray((*PS).rcst0knot,2):*exp(b1:*(*PS).t0) ///
			:+ u:*(asarray((*PS).rcst0knot,1) :- asarray((*PS).rcst0knot,2):*(*PS).t0):*v:*exp(b1:*((*PS).t0)))) :*((*PS).includefirstint)
		}
	}	
	else {
		if (!(*PS).delentry){
			cumhazardscoremin = (u:*v:*(*PS).lowt:^(w):*log((*PS).lowt):*asarray((*PS).rcsfirstknot,2) ///
				:+ u:*(*PS).lowt:^(w):*(-1):*(v:^2):*asarray((*PS).rcsfirstknot,2) ///
				:+ v:*(*PS).lowt:^(w):*u:*(asarray((*PS).rcsfirstknot,1) :- asarray((*PS).rcsfirstknot,2):*log((*PS).lowt))):*((*PS).includefirstint)			
		}
		else {
			cumhazardscoremin = ((u:*v:*log((*PS).lowt):*(*PS).lowt:^(w):*asarray((*PS).rcsfirstknot,2) ///
			:- u:*(v:^2):*(*PS).lowt:^(w):*asarray((*PS).rcsfirstknot,2) ///
			:+ u:*(asarray((*PS).rcsfirstknot,1) :- asarray((*PS).rcsfirstknot,2):*log((*PS).lowt)):*v:*(*PS).lowt:^(w)) ///
			:- (u:*v:*log((*PS).t0):*(*PS).t0:^(w):*asarray((*PS).rcst0knot,2) ///
			:- u:*(v:^2):*(*PS).t0:^(w):*asarray((*PS).rcst0knot,2) ///
			:+ u:*(asarray((*PS).rcst0knot,1) :- asarray((*PS).rcst0knot,2):*log((*PS).t0)):*v:*(*PS).t0:^(w))):*((*PS).includefirstint)
		
		}
	}

// After last knot (analytical)
	if ((*PS).bhtime) {
		cumhazardscoremax = ((u_last:*v_last:*exp(b1_last:*(*PS).t):*(*PS).basedrcs:*(*PS).t ///
		:- u_last:*(v_last:^2):*(*PS).basedrcs:*exp(b1_last:*(*PS).t) ///
		:+ u_last:*((*PS).basercs :- (*PS).basedrcs:*(*PS).t):*v_last:*exp(b1_last:*(*PS).t)) ///
		:- ((u_last:*v_last:*exp(b1_last:*(*PS).hight):*asarray((*PS).rcslastknot,2):*(*PS).hight ///
		:- u_last:*(v_last:^2):*asarray((*PS).rcslastknot,2):*exp(b1_last:*(*PS).hight) ///
		:+ u_last:*(asarray((*PS).rcslastknot,1) :- asarray((*PS).rcslastknot,2):*(*PS).hight):*v_last:*exp(b1_last:*(*PS).hight)))):*((*PS).includethirdint)
	}
	else {
		cumhazardscoremax =((u_last:*v_last:*(*PS).t:^(b1_last:+1):*log((*PS).t):* (*PS).basedrcs ///
			:+ u_last:*(*PS).t:^(b1_last:+1):*(-1):*(v_last:^2):*(*PS).basedrcs:+ ///
			v_last:*(*PS).t:^(b1_last:+1):*u_last:*((*PS).basercs :- (*PS).basedrcs:*log((*PS).t))) ///
			:- (u_last:*v_last:*(*PS).hight:^(b1_last:+1):*log((*PS).hight):*asarray((*PS).rcslastknot,2) ///
			:+ u_last:*(*PS).hight:^(b1_last:+1):*(-1):*(v_last:^2):*asarray((*PS).rcslastknot,2) ///
			:+ v_last:*(*PS).hight:^(b1_last:+1):*u_last:*(asarray((*PS).rcslastknot,1) :- asarray((*PS).rcslastknot,2):*log((*PS).hight)))):*((*PS).includethirdint)	 
		
		}
		
		S = J((*PS).Nobs,0,.)
	
		if(!(*PS).nocons) nsplines = cols(rcsbeta) -1	
		else nsplines =cols(rcsbeta)
		if ((*PS).hascovs) {
			if ((*PS).rsmodel) {
				for (i=1;i<=(*PS).Ncovs;i++) {
					derivcov= ((*PS).d:*(exp((*PS).xb :+ (*PS).rcs):/((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs))) :- cumhazard):* (*PS).covsmat[,i] 									
					S=S,derivcov
				}
			}
			else {
				for (i=1;i<=(*PS).Ncovs;i++) {
					derivcov= ((*PS).d :- (cumhazard)):* (*PS).covsmat[,i]
					S=S,derivcov					
				}
			}
		}

		for (i=1;i<=nsplines;i++) {
			if ((*PS).rsmodel) {
				deriv = (*PS).d:*(exp((*PS).xb :+ (*PS).rcs):/((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs))):*(*PS).basercs[,i]  :- cumhazardscoremin[,i] :- cumhazardscoremax[,i]
				for (j=1;j<=(*PS).Nnodes;j++) {
					deriv = deriv :- (*PS).cumhazardmid[,j] :* (asarray((*PS).nodes, j)[,i])
				}

			}
			else {
				deriv = (*PS).d:*(*PS).basercs[,i]  :- cumhazardscoremin[,i] :- cumhazardscoremax[,i]
				for (j=1;j<=(*PS).Nnodes;j++) {
					deriv = deriv :- (*PS).cumhazardmid[,j] :* (asarray((*PS).nodes, j)[,i])
				}
			}
			S = S,deriv	
		}
		if(!(*PS).nocons) {
			if ((*PS).rsmodel) {
				S = S,((*PS).d:*exp((*PS).xb :+ (*PS).rcs):/((*PS).bhazard :+ exp((*PS).xb :+(*PS).rcs)) :- cumhazard)
			}
			else {
				S = S,((*PS).d :- cumhazard)
			}
		}
		if (todo==1) return	

		

	
////////////////////////////////////////////////////////////////////////////////////////////////		
//HESSIAN


		firstrcs= asarray((*PS).rcsfirstknot,1)[,1]
		firstdrcs =asarray((*PS).rcsfirstknot,2)[,1]
		if ((*PS).delentry) {
			firstrcst0= asarray((*PS).rcst0knot,1)[,1]
			firstdrcst0= asarray((*PS).rcst0knot,2)[,1]
		}
		else {
			firstrcst0= J((*PS).Nobs,1,0)
			firstdrcst0= J((*PS).Nobs,1,0)
		}

// calculate the hessian if there are covariates and splines
	if ((*PS).hascovs) {
		deriv2cov = J((*PS).Nobs,1,0)
		if((*PS).rsmodel) deriv2cov= (*PS).d:*(*PS).covsmat[,1]:*(*PS).covsmat[,1]:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
		// first covariate
		deriv2cov= deriv2cov:- cumhazard:*(*PS).covsmat[,1]:*(*PS).covsmat[,1]
		H= quadcolsum(deriv2cov,1)
		
		//covariate with covariate
		for (i=2;i<=(*PS).Ncovs;i++) {
			ind=1
			t1 = J(1,0,.)
			while (ind<i) {
				deriv2cov = J((*PS).Nobs,1,0)
				if((*PS).rsmodel) deriv2cov= (*PS).d:*(*PS).covsmat[,ind]:*(*PS).covsmat[,i]:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
				deriv2cov = deriv2cov :- cumhazard:*(*PS).covsmat[,ind]:*(*PS).covsmat[,i]
				t1=t1, quadcolsum(deriv2cov,1)
				ind++	
				
			}
			deriv2cov = J((*PS).Nobs,1,0)
			if((*PS).rsmodel) deriv2cov= (*PS).d:*(*PS).covsmat[,i]:*(*PS).covsmat[,i]:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
			deriv2cov = deriv2cov :- cumhazard:*(*PS).covsmat[,i]:*(*PS).covsmat[,i]
			H=H,t1'\t1,quadcolsum(deriv2cov,1)
			
		}
		
		//covariate with splines
		//spline 1 hessian components
		t1=J(1,0,.)
		for (i=1;i<=(*PS).Ncovs;i++) {		
			//second derivatives for spline 1 and covariates 
			if(!(*PS).bhtime){
				deriv2rcscov(covrcs1lowt, firstrcs, firstrcst0, firstdrcs, firstdrcst0, (*PS).covsmat[,i], (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
				deriv2rcscov(covrcs1hight, (*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1], (*PS).covsmat[,i], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
				deriv2rcsrcs(rcssquare1lowt, firstrcs, firstrcst0, firstdrcs, firstdrcst0, firstrcs, firstrcst0, firstdrcs, firstdrcst0, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
				deriv2rcsrcs(rcssquare1hight, (*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1],(*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
			}
			else {		
				deriv2rcscov(covrcs1lowt, firstrcs, firstrcst0, firstdrcs, firstdrcst0, (*PS).covsmat[,i], (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs,0)			
				deriv2rcscov(covrcs1hight, (*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1], (*PS).covsmat[,i], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs,1)			
				deriv2rcsrcs(rcssquare1lowt, firstrcs, firstrcst0, firstdrcs, firstdrcst0, firstrcs, firstrcst0, firstdrcs, firstdrcst0, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
				deriv2rcsrcs(rcssquare1hight, (*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1],(*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
			}
			deriv2cov = J((*PS).Nobs,1,0)		
			if((*PS).rsmodel) {
				for (j=1;j<=(*PS).Nnodes;j++) deriv2cov= (*PS).d:*(asarray((*PS).nodes, j)[,1]):*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2			
			}
			for (j=1;j<=(*PS).Nnodes;j++) deriv2cov = deriv2cov :- (*PS).cumhazardmid[,j]:*(asarray((*PS).nodes, j)[,1])
			deriv2cov =deriv2cov:*(*PS).covsmat[,i]  :- covrcs1lowt :- covrcs1hight
			t1=t1, quadsum(deriv2cov,1)	
			
		}
		deriv2cov = J((*PS).Nobs,1,0)
		if((*PS).rsmodel) {
			for (j=1;j<=(*PS).Nnodes;j++) deriv2cov= (*PS).d:*(asarray((*PS).nodes, j)[,1]):^2:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
		}
		deriv2cov = deriv2cov :- rcssquare1lowt :- rcssquare1hight
		for (j=1;j<=(*PS).Nnodes;j++) deriv2cov= deriv2cov :- (*PS).cumhazardmid[,j]:*(asarray((*PS).nodes, j)[,1]):^2
		H=H,t1'\t1,quadsum(deriv2cov,1)

	// from spline 2 onwards
		for (k=2;k<=nsplines;k++) {
			t1 = J(1,0,.)
			if (((*PS).delentry)) {
				rcst0knotk=asarray((*PS).rcst0knot,1)[,k]
				drcst0knotk=asarray((*PS).rcst0knot,2)[,k]
			}
			else {
				rcst0knotk = J((*PS).Nobs,1,0)
				drcst0knotk = J((*PS).Nobs,1,0)
			}
			for (i=1;i<=(*PS).Ncovs;i++) {
				//second derivatives for spline2+ and covariates 
				if(!(*PS).bhtime){
					deriv2rcscov(covrcslowt, asarray((*PS).rcsfirstknot,1)[,k],rcst0knotk, asarray((*PS).rcsfirstknot,2)[,k],drcst0knotk, (*PS).covsmat[,i], (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)			
					deriv2rcscov(covrcshight, (*PS).basercs[,k], asarray((*PS).rcslastknot,1)[,k], (*PS).basedrcs[,k], asarray((*PS).rcslastknot,2)[,k], (*PS).covsmat[,i], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)			
				}
				else {				
					deriv2rcscov(covrcslowt, asarray((*PS).rcsfirstknot,1)[,k], rcst0knotk, asarray((*PS).rcsfirstknot,2)[,k], drcst0knotk, (*PS).covsmat[,i], (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs,0)			
					deriv2rcscov(covrcshight, (*PS).basercs[,k], asarray((*PS).rcslastknot,1)[,k], (*PS).basedrcs[,k], asarray((*PS).rcslastknot,2)[,k], (*PS).covsmat[,i], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs,1)
				}
				deriv2cov = J((*PS).Nobs,1,0)
				if((*PS).rsmodel) {
					for (j=1;j<=(*PS).Nnodes;j++) deriv2cov= (*PS).d:*(asarray((*PS).nodes, j)[,k]):*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
				}	
				for (j=1;j<=(*PS).Nnodes;j++) deriv2cov = deriv2cov :- (*PS).cumhazardmid[,j]:*(asarray((*PS).nodes, j)[,k])
				deriv2cov = deriv2cov :* (*PS).covsmat[,i] :- covrcslowt :- covrcshight
				t1=t1, quadsum(deriv2cov,1)	
			}						
		// spline - spline hessian components
			ind=1
			while (ind<k) {
				if (((*PS).delentry)) {
					rcst0knotind=asarray((*PS).rcst0knot,1)[,ind]
					drcst0knotind=asarray((*PS).rcst0knot,2)[,ind]
				}
				else {
					rcst0knotind = J((*PS).Nobs,1,0)
					drcst0knotind = J((*PS).Nobs,1,0)
				}				
				if(!(*PS).bhtime){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
					deriv2rcsrcs(rcslowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,ind],  rcst0knotind,  asarray((*PS).rcsfirstknot,2)[,ind],  drcst0knotind, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)					
					deriv2rcsrcs(rcshight,  (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k],  (*PS).basercs[,ind],  asarray((*PS).rcslastknot,1)[,ind],  (*PS).basedrcs[,ind],  asarray((*PS).rcslastknot,2)[,ind], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)				
				}			
				else {	
					deriv2rcsrcs(rcslowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,ind],  rcst0knotind,  asarray((*PS).rcsfirstknot,2)[,ind],  drcst0knotind, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
					deriv2rcsrcs(rcshight,  (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k],  (*PS).basercs[,ind],  asarray((*PS).rcslastknot,1)[,ind],  (*PS).basedrcs[,ind],  asarray((*PS).rcslastknot,2)[,ind], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)				
				}	
				deriv2cov = J((*PS).Nobs,1,0)
				if((*PS).rsmodel) {
					for (j=1;j<=(*PS).Nnodes;j++) deriv2cov= (*PS).d:*(asarray((*PS).nodes, j)[,k]):*(asarray((*PS).nodes, j)[,ind]):*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
				}					
				for (j=1;j<=(*PS).Nnodes;j++) deriv2cov = deriv2cov :- (*PS).cumhazardmid[,j]:*(asarray((*PS).nodes, j)[,k]):*(asarray((*PS).nodes, j)[,ind])			
				deriv2cov = deriv2cov :- rcslowt :- rcshight
				t1=t1,quadsum(deriv2cov,1)	
				ind++					
			}
		// spline squared hessian components
			if(!(*PS).bhtime){
				deriv2rcsrcs(rcssquarelowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
				deriv2rcsrcs(rcssquarehight, (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
			}
			else{			
				 deriv2rcsrcs(rcssquarelowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)				
				deriv2rcsrcs(rcssquarehight, (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
			}	
			deriv2cov = J((*PS).Nobs,1,0)	
			if((*PS).rsmodel) {
				for (j=1;j<=(*PS).Nnodes;j++) deriv2cov= (*PS).d:*asarray((*PS).nodes, j)[,k]:^2:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
			}					
			deriv2cov = deriv2cov :- rcssquarelowt :- rcssquarehight	
			for (j=1;j<=(*PS).Nnodes;j++) deriv2cov = deriv2cov :- (*PS).cumhazardmid[,j]:*asarray((*PS).nodes, j)[,k]:^2
			H=H,t1'\t1,quadsum(deriv2cov,1)	
		}
	
		// add covariate-constant, spline-constant and constant-constant hessian components
		if(!(*PS).nocons) {
			t1 = J(1,0,.)
			// covariate-constant hessian components
			for(i=1;i<=(*PS).Ncovs;i++){
				deriv2cons = J((*PS).Nobs,1,0)
				if((*PS).rsmodel) deriv2cons= (*PS).d:*(*PS).covsmat[,i]:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2				
				deriv2cons = deriv2cons :-cumhazard:*(*PS).covsmat[,i]
				t1=t1,quadcolsum(deriv2cons,1)
			}				
			// spline-constant hessian components
			for(k=1;k<=nsplines;k++) {
				if (((*PS).delentry)) {
					rcst0knotk=asarray((*PS).rcst0knot,1)[,k]
					drcst0knotk=asarray((*PS).rcst0knot,2)[,k]
				}
				else {
					rcst0knotk = J((*PS).Nobs,1,0)
					drcst0knotk = J((*PS).Nobs,1,0)
				}
				deriv2cons = J((*PS).Nobs,1,0)				
				if(!(*PS).bhtime) {
					deriv2rcscons(rcsconslowt, asarray((*PS).rcsfirstknot,1)[,k], rcst0knotk ,asarray((*PS).rcsfirstknot,2)[,k], drcst0knotk, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry, (*PS).includefirstint, (*PS).Nobs, 0)							
					deriv2rcscons(rcsconshight, (*PS).basercs[,k], asarray((*PS).rcslastknot,1)[,k],(*PS).basedrcs[,k], asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry, (*PS).includethirdint, (*PS).Nobs, 1)							
				}
				else{
					deriv2rcscons(rcsconslowt, asarray((*PS).rcsfirstknot,1)[,k], rcst0knotk ,asarray((*PS).rcsfirstknot,2)[,k], drcst0knotk, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry, (*PS).includefirstint, (*PS).Nobs, 0)
					deriv2rcscons(rcsconshight, (*PS).basercs[,k], asarray((*PS).rcslastknot,1)[,k],(*PS).basedrcs[,k], asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry, (*PS).includethirdint, (*PS).Nobs, 1)						
				}
				if((*PS).rsmodel) {
					for (j=1;j<=(*PS).Nnodes;j++) deriv2cons= (*PS).d:* asarray((*PS).nodes, j)[,k]:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
				}					
				for (j=1;j<=(*PS).Nnodes;j++) deriv2cons = deriv2cons :- (*PS).cumhazardmid[,j] :* asarray((*PS).nodes, j)[,k]
				deriv2cons =deriv2cons :-rcsconslowt :-rcsconshight
				t1=t1,quadsum(deriv2cons,1)				
			}			
			//constant-constant hessian component
			deriv2cons = J((*PS).Nobs,1,0)
			if((*PS).rsmodel) deriv2cons= (*PS).d:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2				
			deriv2cons = deriv2cons :- cumhazard
			H=H,t1'\t1,quadsum(deriv2cons,1)				
		}		
	}
	
	// no covariates
	else {	
	// first spline
		if(!(*PS).bhtime){		
			deriv2rcsrcs(rcssquare1lowt, firstrcs, firstrcst0, firstdrcs, firstdrcst0, firstrcs, firstrcst0, firstdrcs, firstdrcst0, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)		
			deriv2rcsrcs(rcssquare1hight, (*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1],(*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)		
		}						
		else{
			deriv2rcsrcs(rcssquare1lowt, firstrcs, firstrcst0, firstdrcs, firstdrcst0, firstrcs, firstrcst0, firstdrcs, firstdrcst0, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)		
			deriv2rcsrcs(rcssquare1hight, (*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1],(*PS).basercs[,1], asarray((*PS).rcslastknot,1)[,1], (*PS).basedrcs[,1], asarray((*PS).rcslastknot,2)[,1], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
		}		
		deriv2 = J((*PS).Nobs,1,0)	
		if((*PS).rsmodel) {
			for (j=1;j<=(*PS).Nnodes;j++) deriv2= (*PS).d :* (asarray((*PS).nodes, j)[,1]):^2:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2 
		}
		deriv2 = deriv2 :- rcssquare1lowt :- rcssquare1hight
		for (j=1;j<=(*PS).Nnodes;j++) deriv2 = deriv2 :- (*PS).cumhazardmid[,j] :* (asarray((*PS).nodes, j)[,1]):^2
		H = quadsum(deriv2)
		
		// #2+ spline	
		for (k=2;k<=nsplines;k++) {
			if (((*PS).delentry)) {
				rcst0knotk=asarray((*PS).rcst0knot,1)[,k]
				drcst0knotk=asarray((*PS).rcst0knot,2)[,k]
			}
			else {
				rcst0knotk = J((*PS).Nobs,1,0)
				drcst0knotk = J((*PS).Nobs,1,0)
			}
			ind = 1
			t1 = J(1,0,.)
			while (ind<k) {
				if (((*PS).delentry)) {
					rcst0knotind=asarray((*PS).rcst0knot,1)[,ind]
					drcst0knotind=asarray((*PS).rcst0knot,2)[,ind]				
				}
				else {
					rcst0knotind = J((*PS).Nobs,1,0)
					drcst0knotind = J((*PS).Nobs,1,0)
				}
				if(!(*PS).bhtime){
					deriv2rcsrcs(rcslowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,ind],  rcst0knotind,  asarray((*PS).rcsfirstknot,2)[,ind],  drcst0knotind, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
					deriv2rcsrcs(rcshight,  (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k],  (*PS).basercs[,ind],  asarray((*PS).rcslastknot,1)[,ind],  (*PS).basedrcs[,ind],  asarray((*PS).rcslastknot,2)[,ind], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)				
				}
				else {
					deriv2rcsrcs(rcslowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,ind],  rcst0knotind,  asarray((*PS).rcsfirstknot,2)[,ind],  drcst0knotind, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
					deriv2rcsrcs(rcshight,  (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k],  (*PS).basercs[,ind],  asarray((*PS).rcslastknot,1)[,ind],  (*PS).basedrcs[,ind],  asarray((*PS).rcslastknot,2)[,ind], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)							
				}
				deriv2 = J((*PS).Nobs,1,0)
				if((*PS).rsmodel) {
					for (j=1;j<=(*PS).Nnodes;j++) deriv2= (*PS).d :* (asarray((*PS).nodes, j)[,k] :* asarray((*PS).nodes, j)[,ind]):*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2 
				}
				deriv2 = deriv2 :- rcslowt :-rcshight
				for (j=1;j<=(*PS).Nnodes;j++) deriv2 = deriv2 :- (*PS).cumhazardmid[,j] :* (asarray((*PS).nodes, j)[,k] :* asarray((*PS).nodes, j)[,ind])
				t1 = t1,quadsum(deriv2)
				ind++		
			}
			if(!(*PS).bhtime){
				deriv2rcsrcs(rcssquarelowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
				deriv2rcsrcs(rcssquarehight, (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
			}
			else{	
				deriv2rcsrcs(rcssquarelowt,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk,  asarray((*PS).rcsfirstknot,1)[,k],  rcst0knotk,  asarray((*PS).rcsfirstknot,2)[,k],  drcst0knotk, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry,(*PS).includefirstint, (*PS).Nobs, 0)
				deriv2rcsrcs(rcssquarehight, (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).basercs[,k],  asarray((*PS).rcslastknot,1)[,k],  (*PS).basedrcs[,k],  asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry,(*PS).includethirdint, (*PS).Nobs, 1)
			}				
			deriv2 = J((*PS).Nobs,1,0)
			if((*PS).rsmodel) {
				for (j=1;j<=(*PS).Nnodes;j++) deriv2= (*PS).d :* (asarray((*PS).nodes, j)[,k]):^2:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2 
			}
			deriv2 = deriv2 :- rcssquarehight :-rcssquarelowt
			for (j=1;j<=(*PS).Nnodes;j++) deriv2 = deriv2 :- (*PS).cumhazardmid[,j] :* (asarray((*PS).nodes, j)[,k]):^2
			H = H,t1'\t1,quadsum(deriv2)
		}

		// add covariate-constant, spline-constant and constant-constant hessian components
		if(!(*PS).nocons) {
			t1 = J(1,0,.)	
			// spline-constant hessian components
			for(k=1;k<=nsplines;k++) {
				if (((*PS).delentry)) {
					rcst0knotk=asarray((*PS).rcst0knot,1)[,k]
					drcst0knotk=asarray((*PS).rcst0knot,2)[,k]
				}
				else {
					rcst0knotk = J((*PS).Nobs,1,0)
					drcst0knotk = J((*PS).Nobs,1,0)
				}
				if(!(*PS).bhtime) {
					deriv2rcscons(rcsconslowt, asarray((*PS).rcsfirstknot,1)[,k], rcst0knotk ,asarray((*PS).rcsfirstknot,2)[,k], drcst0knotk, (*PS).lowt, (*PS).t0, u, v, w, (*PS).bhtime, (*PS).delentry, (*PS).includefirstint, (*PS).Nobs, 0)							
					deriv2rcscons(rcsconshight, (*PS).basercs[,k], asarray((*PS).rcslastknot,1)[,k],(*PS).basedrcs[,k], asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, w_last, (*PS).bhtime, (*PS).delentry, (*PS).includethirdint, (*PS).Nobs, 1)							
				}
				else{
					deriv2rcscons(rcsconslowt, asarray((*PS).rcsfirstknot,1)[,k], rcst0knotk ,asarray((*PS).rcsfirstknot,2)[,k], drcst0knotk, (*PS).lowt, (*PS).t0, u, v, b1, (*PS).bhtime, (*PS).delentry, (*PS).includefirstint, (*PS).Nobs, 0)							
					deriv2rcscons(rcsconshight, (*PS).basercs[,k], asarray((*PS).rcslastknot,1)[,k],(*PS).basedrcs[,k], asarray((*PS).rcslastknot,2)[,k], (*PS).t, (*PS).hight, u_last, v_last, b1_last, (*PS).bhtime, (*PS).delentry, (*PS).includethirdint, (*PS).Nobs, 1)							
				}
				deriv2cons = J((*PS).Nobs,1,0)				
				if((*PS).rsmodel) {
					for (j=1;j<=(*PS).Nnodes;j++) deriv2cons= (*PS).d:* asarray((*PS).nodes, j)[,k]:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
				}
				for (j=1;j<=(*PS).Nnodes;j++) deriv2cons = deriv2cons :- (*PS).cumhazardmid[,j] :* asarray((*PS).nodes, j)[,k]
				deriv2cons =deriv2cons :-rcsconslowt :-rcsconshight
				t1=t1,quadsum(deriv2cons,1)				
			}				
			//constant-constant hessian component
			deriv2cons = J((*PS).Nobs,1,0)
			if((*PS).rsmodel) deriv2cons= (*PS).d:*(*PS).bhazard:*exp((*PS).xb :+ (*PS).rcs) :/ ((*PS).bhazard :+ exp((*PS).xb :+ (*PS).rcs)):^2
			deriv2cons = deriv2cons :- cumhazard
			H=H,t1'\t1,quadsum(deriv2cons,1)				
		}			
	}
	if (todo==2) return		
}
	
	
////////////////////////////////////////////////////////////////////////////////////////////////		
// Functions for Hessian quadrature	
// Spline-covariate parts				
function deriv2rcscov(created, splinet1, splinet2, dsplinet1, dsplinet2, matrixname, time1, time2, u, v, w, bhtime, delentry, interval, nobs, hight)
{
	if(hight==0){
		if(!(bhtime)) { 
			if (!(delentry)) {
				created = ((u:*v:*(time1):^w:*matrixname:*(splinet1:-v:*dsplinet1))):*(interval)
			}
			else {
				created = ((u:*v:*(time1):^w:*matrixname:*(splinet1 :- v:*dsplinet1)) :-  ///
					(u:*v:*(time2):^w:*matrixname:*(splinet2 :- v:*dsplinet2))):*(interval)
			}
		}
		else {
			if(!(delentry)) {
				created = ((u:*v:*matrixname:*(exp(w:*time1):*(splinet1 :- v:*dsplinet1) :+ dsplinet1:*(v :+ time1) :-splinet1)))				
			}
			else {
				created = ((u:*v:*matrixname:*(exp(w:*time1):*(splinet1 :- v:* dsplinet1) :+ (v :+ time1):*dsplinet1 :- splinet1)) :- ///
				((u:*v:*matrixname:*(exp(w:*time2):*(splinet2 :- v:* dsplinet2) :+ (v :+ time2):*dsplinet2 :- splinet2)))):*(interval)
			}
		}
	}
	if (hight==1) {
		if (!(bhtime)) {
			created = ((u:*v:*(time1):^w:*matrixname:*(splinet1 :- v:*dsplinet1)) :-  ///
					(u:*v:*(time2):^w:*matrixname:*(splinet2 :- v:*dsplinet2))):*(interval)
		}
		else {
				created = ((u:*v:*matrixname:*(exp(w:*time1):*(splinet1 :- v:* dsplinet1) :+ (v :+ time1):*dsplinet1 :- splinet1)) :- ///
				((u:*v:*matrixname:*(exp(w:*time2):*(splinet2 :- v:* dsplinet2) :+ (v :+ time2):*dsplinet2 :- splinet2)))):*(interval)
		}
	}
}

// Spline-spline parts
function deriv2rcsrcs(created, spline1, spline1t2, dspline1, dspline1t2, spline2, spline2t2, dspline2, dspline2t2, time1, time2, u, v, w, bhtime, delentry, interval, nobs, hight)
{
		if(hight==0){
			if(!(bhtime)) { 
				if (!(delentry)) {
					created = (u:*v:*time1:^w:*(dspline1:*(2:*v:^2:*dspline2 :- spline2:*v) :+ spline1:*(spline2 :- v:*dspline2))):*(interval)	
				}
				else {
					created = ((u:*v:*time1:^w:*(dspline1:*(2:*v:^2:*dspline2 :- spline2:*v) :+ spline1:*(spline2 :- v:*dspline2))):- ///
								(u:*v:*time2:^w:*(dspline1t2:*(2:*v:^2:*dspline2t2 :- spline2t2:*v) :+ spline1t2:*(spline2t2 :- v:*dspline2t2)))):*(interval)
				}
			}
			else {
				if(!(delentry)) {
					created = (u:*v:*exp(w:*time1):*(dspline1:*(2:*v:^2:*dspline2 :- spline2:*v) :+ spline1:*(spline2:-v:*dspline2)) :- ///
								u:*v:*(dspline1:*dspline2:*(2:*v:^2:+2:*v:*time1:+(time1):^2) :- dspline1:*spline2:*(v:+time1) :- dspline2:*spline1:*(v:+time1) + spline1:*spline2)):*(interval)
				}
				else {
					created = ((u:*v:*exp(w:*time1):*(dspline1:*(2:*v:^2:*dspline2 :- spline2:*v) :+ spline1:*(spline2 :- v:*dspline2))) :- ///
						(u:*v:*exp(w:*time2):*(dspline1t2:*(2:*v:^2:*dspline2t2 :- spline2t2:*v) :+ spline1t2:*(spline2t2 :- v:*dspline2t2)))) :* (interval)
				}
			}
		}
		if(hight==1) {
			if (!(bhtime)) {
				created = ((u:*v:*time1:^w:* (dspline1:*(2:*v:^2:*dspline2 :- spline2:*v) :+ spline1:*(spline2 :- v:*dspline2))) :- ///
				(u:*v:*time2:^w:* (dspline1t2:*(2:*v:^2:*dspline2t2 :- spline2t2:*v) :+ spline1t2:*(spline2t2 :- v:*dspline2t2)))):*(interval)	
			}
			else {
				created = ((u:*v:*exp(w:*time1):*(dspline1:*(2:*v:^2:*dspline2 :- spline2:*v) :+ spline1 :*(spline2 :- v:*dspline2))) :- ///
					(u:*v:*exp(w:*time2):*(dspline1t2:*(2:*v:^2:*dspline2t2 :- spline2t2:*v) :+ spline1t2 :*(spline2t2 :- v:*dspline2t2)))):*(interval)
			}
		}
}



//spline-constant
function deriv2rcscons(created, spline1, spline1t2, dspline1, dspline1t2, time1, time2, u, v, w, bhtime, delentry, interval, nobs, hight)
{
		if(hight==0){
			if(!(bhtime)) { 
				if (!(delentry)) {
					created = (u:*v:*(time1):^w:*(spline1 :- v:*dspline1 )):*(interval)	
				}
				else {
					created = ((u:*v:*(time1):^w:*(spline1 :- v:*dspline1 )):- (u:*v:*(time2):^w:*(spline1t2 :- v:*dspline1t2 ))):*(interval)					
				}
			}
			else {
				if(!(delentry)) {
					created = (u:*v:*(exp(w:*(time1)):*(spline1 :- v:*dspline1) :+ dspline1:*(v :+ time1) :- spline1)):*(interval)
				}
				else {
					//created = ((u:*v:*(exp(w:*(time1)):*(spline1 :- v:*dspline1) :+ dspline1:*(v :+ time1) :- spline1)):- ///
						//(u:*v:*(exp(w:*(time2)):*(spline1t2 :- v:*dspline1t2) :+ dspline1t2:*(v :+ time2) :- spline1t2))):*(interval)
					created = (u:*v:*((exp(w:*time1):*(spline1 :- v:*dspline1)) :- exp(w:*time2):*(spline1t2 :- v:*dspline1t2))):*(interval)
				}
			}
		}
		if (hight==1) {
			if (!(bhtime)) {
				created  = ((u:*v:*(time1):^w:*(spline1 :- v:*dspline1 )) :- (u:*v:*(time2):^w:*(spline1t2 :- v:*dspline1t2 ))):*(interval)
			}
			else {
				created = ((u:*v:*(exp(w:*(time1)):*(spline1 :-v:*dspline1)))  :-(u:*v:*exp(w:*time2):*(spline1t2 :- v:*dspline1t2))) :* (interval)
			}
		}
}


// Gaussian Quadrature 
void gq(string scalar weightsname, string scalar nodesname)
{
        n =  strtoreal(st_local("n"))
        inttype = st_local("inttype")
        i = range(1,n,1)'
        i1 = range(1,n-1,1)'
        alpha = strtoreal(st_local("alpha"))
        beta = strtoreal(st_local("beta"))
                
        if(inttype == "legendre") {
                muzero = 2
                a = J(1,n,0)
                b = i1:/sqrt(4 :* i1:^2 :- 1)
        }
        else if(inttype == "cheb1") {
                muzero = pi()
                a = J(1,n,0)
                b = J(1,n-1,0.5)
                b[1] = sqrt(0.5)
    }
        else if(inttype == "cheb2") {
                muzero = pi()/2
                a = J(1,n,0)
                b = J(1,n-1,0.5)
        }
        else if(inttype == "hermite") {
                muzero = sqrt(pi())
                a = J(1,n,0)
                b = sqrt(i1:/2)
        }
        else if(inttype == "jacobi") {
                ab = alpha + beta
                muzero = 2:^(ab :+ 1) :* gamma(alpha :+ 1) * gamma(beta :+ 1):/gamma(ab :+ 2)
                a = i
                a[1] = (beta - alpha):/(ab :+ 2)
                i2 = range(2,n,1)'
                abi = ab :+ (2 :* i2)
                a[i2] = (beta:^2 :- alpha^2):/(abi :- 2):/abi
                b = i1
        b[1] = sqrt(4 * (alpha + 1) * (beta + 1):/(ab :+ 2):^2:/(ab :+ 3))
        i2 = i1[2..n-1]
        abi = ab :+ 2 :* i2
        b[i2] = sqrt(4 :* i2 :* (i2 :+ alpha) :* (i2 :+ beta) :* (i2 :+ ab):/(abi:^2 :- 1):/abi:^2)
        }
        else if(inttype == "laguerre") {
                a = 2 :* i :- 1 :+ alpha
                b = sqrt(i1 :* (i1 :+ alpha))
                muzero = gamma(alpha :+ 1)
    }

        A= diag(a)
        for(j=1;j<=n-1;j++){
                A[j,j+1] = b[j]
                A[j+1,j] = b[j]
        }
        symeigensystem(A,vec,nodes)
        weights = (vec[1,]:^2:*muzero)'
        weights = weights[order(nodes',1)]
        nodes = nodes'[order(nodes',1)']
        st_matrix(weightsname,weights)
        st_matrix(nodesname,nodes)
}
                
end

