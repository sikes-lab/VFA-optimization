// This macro is to be used in conjunction with the MATLAB scripts 
// 'Process_Cyan_WhiteBkgrd.m' or 'Process_Cyan_BlackBkgrd.m'

First_window = getTitle;

run("RGB to CMYK");
run("Stack to Images");
wait(2);

//Start closing unwanted windows
selectWindow(First_window);
close();
selectWindow("M");
close();
selectWindow("Y");
close();
selectWindow("K");
close();
//end close

selectWindow("C");
run("Measure");
close();

//wait(1);

//selectWindow("Results");
//run("Close");
