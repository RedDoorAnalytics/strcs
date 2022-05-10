
capture erase lstrcs.mlib
mata: mata set matastrict off
qui {
        do "./strcs/strcs_matacode.mata"
        
        mata: mata mlib create lstrcs, dir(.)
        mata: mata mlib add    lstrcs *(), dir(.)
        mata: mata d *()  
        mata mata clear
        mata mata mlib index
}
