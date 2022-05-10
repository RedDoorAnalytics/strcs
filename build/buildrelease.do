//build new version of strcs
// --> run whole do file

local drive /Users/Michael/Documents/reddooranalytics/products/strcs
cd `drive'

//=======================================================================================================================//
	
//build for SSC -> current version up is 1.85
	
local sscversion 1_85
cap mkdir ./ssc/version_`sscversion'
local fdir `drive'/ssc/version_`sscversion'/

//=======================================================================================================================//

//pkg files
copy ./build/README.txt `fdir', replace
	
//=======================================================================================================================//

//strcs

copy ./strcs/strcs.ado `fdir', replace
copy ./strcs/strcs_pred.ado `fdir', replace
copy ./lstrcs.mlib `fdir', replace

//help files
copy ./strcs/strcs.sthlp `fdir', replace
copy ./strcs/strcs_postestimation.sthlp `fdir', replace
	
