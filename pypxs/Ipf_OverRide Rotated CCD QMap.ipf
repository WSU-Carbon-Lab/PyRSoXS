#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method.

////////////////////////////////////////////////////////////////////////////////////////////
//
// Written by Devin Grabner, Washington State University, Nov 2023
//					********* LAST UPDATED: 24 Jan 2024  *************
//
//	If you have any further questions please contact Devin Grabner (devin.grabner@wsu.edu) or 
//		Brian Collins (brian.collins@wsu.edu) for more information.
//
// Function NI1A_Create2DQWave(DataWave) is from NIKA "NI1_ConvProc" 
// Use this to overload the NIKA function for calculating the QMap for the CCD
//		image and applying an appropriate Solid Angle Correction. This function allows for a rotated
//		CCD Geometry where the sample doesn't have to be located at the axis of rotation of the CCD.
//
// The additional parameters needed are not currently in NIKA
//
// d1 = Distance from axis of rotation to sample (mm)
//				A postive number is assuming the sample is in front of the CCD axis of rotation

// d2 = Sample - Detector distance (mm)
//			This could be pulled from the Sample Detector Distance in the 2D to 1D Data
//			Conversion Panel but I didn't set it to pull from there until 'd1' is also an
//			option with it.
//
////////////////////////////////////////////////////////////////////////////////////////////


Override Function NI1A_Create2DQWave(DataWave)
	Wave DataWave
	string OldDf=GetDataFolder(1)
	
	Wave/Z Param_Values = root:Param_Values
	
	if(WaveExists(Param_Values) == 0)
		print "Making Needed Paramater Waves with Default Values"
		setDataFolder root
		Make/T/N=2 Param_Names = {"Distance from axis of rotation to sample (mm)", "Sample - Detector distance (mm)"}
		Make/S/N=2 Param_Values = {26, 23}
		
		Wave/Z Param_Values = root:Param_Values
	endif
	
	Variable d1 = Param_Values[0] // Distance from axis of rotation to sample (mm)
	Variable d2 = Param_Values[1] // Sample - Detector distance (mm)

	SVAR header = root:Packages:Nika_RSoXS:NRS_Vars:header
	variable phi = str2num(Stringbykey("CCD Theta",header,"=",";")) // Angle the CCD is rotated to (Input: degrees)
	
	// Because of the tan(x) and arctan(x) function it has a hard time handeling exactly 0 and 90 degrees.
	// so this statment moves it just slightly of so the calculation works.
	if (phi == 0)
		phi = (1e-11)
	elseif (phi == 90)
		phi = 89.9999999
	endif
	
	phi = phi * (pi/180) // Convert to radians
	
	//IN2G_PrintDebugStatement(IrenaDebugLevel, 5,"")
	setDataFolder root:Packages:Convert2Dto1D
	
	Wave/Z Q2DWave
	string NoteStr
	NoteStr = note(DataWave)
	NVAR Wavelength = root:Packages:Convert2Dto1D:Wavelength							//in A
	NVAR PixelSizeX = root:Packages:Convert2Dto1D:PixelSizeX							//in millimeters
	NVAR PixelSizeY = root:Packages:Convert2Dto1D:PixelSizeY							//in millimeters
	NVAR beamCenterX=root:Packages:Convert2Dto1D:beamCenterX
	NVAR beamCenterY=root:Packages:Convert2Dto1D:beamCenterY
	
	Variable C = BeamCenterY*PixelSizeY	//Distance from the beam center in normal incidence to bottom of the CCD (mm)

	print "Creating 2D Q wave"
	//Create wave for q distribution
	MatrixOp/O/NTHR=0 Q2DWave=DataWave
	Redimension/S Q2DWave

	Duplicate/FREE Q2DWave, a, b
	
	multithread a = (d1+d2)*sec(phi)-d1-(((q+0.5)*PixelSizeY)+(d1+d2)*tan(phi)-C)*sin(phi)		
	multithread b = ((((q+0.5)*PixelSizeY)+(d1+d2)*tan(phi)-C)*cos(phi))^2 + ((p - BeamCenterX + 0.5)*PixelSizeX)^2

	MatrixOP/O Q2DWave = (2*pi*sqrt(2)/Wavelength)*sqrt(1-a/sqrt(sq(a)+b))
	
	NoteStr = ReplaceStringByKey("BeamCenterX", NoteStr, num2str(BeamCenterX), "=", ";")
	NoteStr = ReplaceStringByKey("BeamCenterY", NoteStr, num2str(BeamCenterY), "=", ";")
	NoteStr = ReplaceStringByKey("PixelSizeX", NoteStr, num2str(PixelSizeX), "=", ";")
	NoteStr = ReplaceStringByKey("PixelSizeY", NoteStr, num2str(PixelSizeY), "=", ";")
	NoteStr = ReplaceStringByKey("Wavelength", NoteStr, num2str(Wavelength), "=", ";")
	note/K Q2DWave, NoteStr
	setDataFolder OldDf
end