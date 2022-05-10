local drive /Users/Michael/Documents/reddooranalytics/products/strcs
cd `drive'

adopath ++ "./strcs"
pr drop _all

// do ./build/buildmlib.do

webuse brcancer,clear
stset rectime, f(censrec)

strcs hormon, df(3)
