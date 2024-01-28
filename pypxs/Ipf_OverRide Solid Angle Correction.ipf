#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

////////////////////////////////////////////////////////////////////////////////////////////
// Because of the way that NIKA does prefactor corrections "NI1_ConvProc.ipf" line #577
//		The solid angle correction has to be done on the 'Calibrated2DDataSet' in
//		the Function NI1A_CorrectDataPerUserReq(orientation) which is copied and 'Override' below
//
//	Function/WAVE SolidAngleCorrection(DataWave) was written by
//			Devin Grabner, Washington State University, 24 Jan 2024
//					********* LAST UPDATED: N/A  *************
//
//	If you have any further questions please contact Devin Grabner (devin.grabner@wsu.edu) or 
//		Brian Collins (brian.collins@wsu.edu) for more information.
//
// Needs to be used in conjunction with the "Overloaded Rotated CCD QMap"
//		(Override Function NI1A_Create2DQWave(DataWave)) from NIKA "NI1_ConvProc" 
// 
//	Use this to overload the NIKA function for calculating an appropriate Solid Angle Correction.
//		This function allows for a rotated CCD Geometry where the sample doesn't have to be
//		located at the axis of rotation of the CCD.
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

Override Function NI1A_CorrectDataPerUserReq(orientation)
	string orientation
	//IN2G_PrintDebugStatement(IrenaDebugLevel, 5,"")
	string OldDf = GetDataFolder(1)
	setDataFolder root:Packages:Convert2Dto1D

	SVAR CalibrationFormula=root:Packages:Convert2Dto1D:CalibrationFormula
	NVAR UseSampleThickness= root:Packages:Convert2Dto1D:UseSampleThickness
	NVAR UseSampleTransmission= root:Packages:Convert2Dto1D:UseSampleTransmission
	NVAR UseCorrectionFactor= root:Packages:Convert2Dto1D:UseCorrectionFactor
	NVAR UseSolidAngle=root:Packages:Convert2Dto1D:UseSolidAngle
	NVAR UseMask= root:Packages:Convert2Dto1D:UseMask
	NVAR UseDarkField= root:Packages:Convert2Dto1D:UseDarkField
	NVAR UseEmptyField= root:Packages:Convert2Dto1D:UseEmptyField
	NVAR UseSubtractFixedOffset= root:Packages:Convert2Dto1D:UseSubtractFixedOffset
	NVAR UseSampleMeasTime= root:Packages:Convert2Dto1D:UseSampleMeasTime
	NVAR UseEmptyMeasTime= root:Packages:Convert2Dto1D:UseEmptyMeasTime
	NVAR UseDarkMeasTime= root:Packages:Convert2Dto1D:UseDarkMeasTime
	NVAR UsePixelSensitivity= root:Packages:Convert2Dto1D:UsePixelSensitivity
	NVAR UseI0ToCalibrate = root:Packages:Convert2Dto1D:UseI0ToCalibrate
	NVAR UseMonitorForEF = root:Packages:Convert2Dto1D:UseMonitorForEF
	
	NVAR PixelSizeX=root:Packages:Convert2Dto1D:PixelSizeX
	NVAR PixelSizeY=root:Packages:Convert2Dto1D:PixelSizeY
	NVAR SampleToCCDDistance=root:Packages:Convert2Dto1D:SampleToCCDDistance
	
	NVAR CorrectionFactor = root:Packages:Convert2Dto1D:CorrectionFactor
	NVAR SampleI0 = root:Packages:Convert2Dto1D:SampleI0
	NVAR EmptyI0 = root:Packages:Convert2Dto1D:EmptyI0
	NVAR SampleThickness = root:Packages:Convert2Dto1D:SampleThickness
	NVAR SampleTransmission = root:Packages:Convert2Dto1D:SampleTransmission
	NVAR SampleMeasurementTime = root:Packages:Convert2Dto1D:SampleMeasurementTime
	NVAR BackgroundMeasTime = root:Packages:Convert2Dto1D:BackgroundMeasTime
	NVAR EmptyMeasurementTime = root:Packages:Convert2Dto1D:EmptyMeasurementTime
	NVAR SubtractFixedOffset = root:Packages:Convert2Dto1D:SubtractFixedOffset
	NVAR CorrectSelfAbsorption = root:Packages:Convert2Dto1D:CorrectSelfAbsorption
	
	NVAR DoGeometryCorrection=root:Packages:Convert2Dto1D:DoGeometryCorrection

	NVAR Use2DPolarizationCor = root:Packages:Convert2Dto1D:Use2DPolarizationCor
	NVAR DoPolarizationCorrection = root:Packages:Convert2Dto1D:DoPolarizationCorrection

	NVAR UseCalib2DData=root:Packages:Convert2Dto1D:UseCalib2DData
	
	Wave DataWave=root:Packages:Convert2Dto1D:CCDImageToConvert
	Wave/Z EmptyRunWave=root:Packages:Convert2Dto1D:EmptyData
	Wave/Z DarkCurrentWave=root:Packages:Convert2Dto1D:DarkFieldData
	Wave/Z MaskWave=root:Packages:Convert2Dto1D:M_ROIMask
	Wave/Z Pix2DSensitivity=root:Packages:Convert2Dto1D:Pixel2DSensitivity
	//little checking here...
	if(UseMask)
		if(!WaveExists(MaskWave) || DimSize(MaskWave, 0)!=DimSize(DataWave, 0) || DimSize(MaskWave, 1)!=DimSize(DataWave, 1))
			abort "Mask problem - either does not exist or has differnet dimensions that data "
		endif
	endif
	if(UseDarkField&&!UseCalib2DData)
		if(!WaveExists(DarkCurrentWave) || DimSize(DarkCurrentWave, 0)!=DimSize(DataWave, 0) || DimSize(DarkCurrentWave, 1)!=DimSize(DataWave, 1))
			abort "Dark field problem - either does not exist or has differnet dimensions that data "
		endif
	endif

	Duplicate/O DataWave,  Calibrated2DDataSet, CalibrationPrefactor
	
	Wave Calibrated2DDataSet=root:Packages:Convert2Dto1D:Calibrated2DDataSet
	redimension/S Calibrated2DDataSet
	string OldNote=note(Calibrated2DDataSet)

	variable tempVal
	MatrixOP/O CalibrationPrefactor = CalibrationPrefactor * 0
	MatrixOP/O CalibrationPrefactor = CalibrationPrefactor + 1 // Sets the whole matrix equal to 1
	
	if(!UseCalib2DData)
		if(UseCorrectionFactor)
			//MatrixOP/O CalibrationPrefactor = CalibrationPrefactor * CorrectionFactor
		endif
		if(UseSolidAngle)
			Duplicate/O CalibrationPrefactor, SolidAngle
			WAVE SolidAngle = SolidAngleCorrection(CalibrationPrefactor)
			MatrixOP/O CalibrationPrefactor = CalibrationPrefactor*SolidAngle
		endif
		if(UseI0ToCalibrate)
			MatrixOP/O CalibrationPrefactor = CalibrationPrefactor / SampleI0
		endif
		if(UseSampleThickness)
			MatrixOP/O CalibrationPrefactor = CalibrationPrefactor / (SampleThickness/10)		//NOTE: changed in ver 1.75 (1/2017), previously was not converted to cm.
			//this is potentially breaking calibration of prior experiments. Need User warning! 
		endif
	
		Duplicate/O DataWave, tempDataWv, tempEmptyField
		redimension/S tempDataWv, tempEmptyField
		
		if(UsePixelSensitivity)
			MatrixOP/O  tempDataWv=tempDataWv/Pix2DSensitivity
		endif
		if(UseDarkField)
			if(UseSampleMeasTime && UseDarkMeasTime)
				if(UsePixelSensitivity)
					tempVal = SampleMeasurementTime/BackgroundMeasTime
					MatrixOP/O  tempDataWv = tempDataWv - (tempVal*DarkCurrentWave/Pix2DSensitivity)
				else
					tempVal = SampleMeasurementTime/BackgroundMeasTime
					MatrixOP/O  tempDataWv = tempDataWv - (tempVal*DarkCurrentWave)
				endif
			else
				if(UsePixelSensitivity)
					MatrixOP/O  tempDataWv = tempDataWv - (DarkCurrentWave/Pix2DSensitivity)
				else
					MatrixOP/O  tempDataWv = tempDataWv - DarkCurrentWave
				endif
			endif
		endif
		if(UseSubtractFixedOffset)
			MatrixOP/O  tempDataWv = tempDataWv - SubtractFixedOffset
		endif
		if(UseSampleTransmission)
			//this is normal correcting by one transmission. 
			MatrixOP/O  tempDataWv=tempDataWv/SampleTransmission
			if(CorrectSelfAbsorption && SampleTransmission<1)
				variable MuCalc=-1*ln(SampleTransmission)/SampleThickness
				variable muD = MuCalc*SampleThickness
				Wave Theta2DWave = root:Packages:Convert2Dto1D:Theta2DWave		//this is actually Theta in radians. 
				if(DimSize(Theta2DWave, 0 )!=DimSize(tempDataWv, 0) || DimSize(Theta2DWave, 1)!=DimSize(tempDataWv, 1) )
					NI1A_Create2DQWave(tempDataWv)				//creates 2-D Q wave does not need to be run always...
					NI1A_Create2DAngleWave(tempDataWv)			//creates 2-D Azimuth Angle wave does not need to be run always...
				endif 	
				MatrixOP/free  SelfAbsorption2D=tempDataWv
				//next is formula 29, chapter 3.4.7 Brain Pauw paper
				MatrixOp/Free  MuDdivCos2TH=MuD*rec(cos(2*Theta2DWave))
				MatrixOp/Free  OneOverBottomPart = rec( -1*MuDdivCos2TH + MuD)
				variable expmud=exp(muD)
				variable expNmud=exp(-1*MuD)
				MatrixOP/O  SelfAbsorption2D=expmud * (exp(-MuDdivCos2TH) - expNmud) * OneOverBottomPart
				//replace nans around center... 
				MatrixOP/O  SelfAbsorption2D=replaceNaNs(SelfAbsorption2D,1)
				//and now correct... 
				MatrixOP/O  tempDataWv=tempDataWv / SelfAbsorption2D
				if(IrenaDebugLevel>1)
					variable MaxCorrection
					wavestats SelfAbsorption2D
					MaxCorrection = 1/wavemin(SelfAbsorption2D)
					print "Sample self absorption correction max is : "+num2str(MaxCorrection) 
				endif
			else
				//print "Could not do corection for self absorption, wrong parameters" 
			endif
		endif
		tempEmptyField=0
		variable ScalingConstEF=1
	
		if(UseEmptyField)
			tempEmptyField = EmptyRunWave
			if(UsePixelSensitivity)
				MatrixOP/O  tempEmptyField = tempEmptyField/Pix2DSensitivity
			endif
			if(UseSubtractFixedOffset)
				MatrixOP/O  tempEmptyField = tempEmptyField - SubtractFixedOffset
			endif
		
			if(UseMonitorForEF)
				ScalingConstEF=SampleI0/EmptyI0
			elseif(UseEmptyMeasTime && UseSampleMeasTime)
				ScalingConstEF=SampleMeasurementTime/EmptyMeasurementTime
			endif
	
			if(UseDarkField)
				if(UseSampleMeasTime && UseEmptyMeasTime)
					if(UsePixelSensitivity)
						tempVal = EmptyMeasurementTime/BackgroundMeasTime
						MatrixOP/O  tempEmptyField=tempEmptyField - (tempVal*(DarkCurrentWave/Pix2DSensitivity))
					else
						tempVal = EmptyMeasurementTime/BackgroundMeasTime
						MatrixOP/O  tempEmptyField=tempEmptyField - (tempVal*DarkCurrentWave)
					endif
				else
					if(UsePixelSensitivity)
						MatrixOP/O  tempEmptyField=tempEmptyField - (DarkCurrentWave/Pix2DSensitivity)
					else
						MatrixOP/O  tempEmptyField=tempEmptyField - DarkCurrentWave
					endif
				endif
			endif
	
		endif
	
		MatrixOP/O  Calibrated2DDataSet = CalibrationPrefactor * (tempDataWv - ScalingConstEF * tempEmptyField)
		
		if(DoGeometryCorrection)  		//geometry correction (= cos(angle)^3) for solid angle projection, added 6/24/2006 to do in 2D data, not in 1D as done (incorrectly also) before using Dales routine.
			NI1A_GenerateGeometryCorr2DWave()
			Wave GeometryCorrection
			MatrixOp/O  Calibrated2DDataSet = Calibrated2DDataSet / GeometryCorrection
		endif
		
//		if(DoPolarizationCorrection)		//added 8/31/09 to enable 2D corection for polarization
//			NI1A_Generate2DPolCorrWv()
//			Wave polar2DWave
//			MatrixOp/O  Calibrated2DDataSet = Calibrated2DDataSet / polar2DWave 		//changed to "/" on October 12 2009 since due to use MatrixOp in new formula the calculate values are less than 1 and this is now correct. 
//		endif
		
		//Add to note:
		//need to add also geometry parameters
		NVAR BeamCenterX=root:Packages:Convert2Dto1D:BeamCenterX
		NVAR BeamCenterY=root:Packages:Convert2Dto1D:BeamCenterY
		NVAR BeamSizeX=root:Packages:Convert2Dto1D:BeamSizeX
		NVAR BeamSizeY=root:Packages:Convert2Dto1D:BeamSizeY
		NVAR HorizontalTilt=root:Packages:Convert2Dto1D:HorizontalTilt
		NVAR XrayEnergy=root:Packages:Convert2Dto1D:XrayEnergy
		NVAR VerticalTilt=root:Packages:Convert2Dto1D:VerticalTilt
		NVAR PixelSizeX=	root:Packages:Convert2Dto1D:PixelSizeX
		NVAR PixelSizeY=	root:Packages:Convert2Dto1D:PixelSizeY
		NVAR SampleToCCDDistance=	root:Packages:Convert2Dto1D:SampleToCCDDistance
		NVAR Wavelength=	root:Packages:Convert2Dto1D:Wavelength
		OldNote+= "Nika_SampleToDetectorDistacne="+num2str(SampleToCCDDistance)+";"
		OldNote+= "Nika_Wavelength="+num2str(Wavelength)+";"
		OldNote+= "Nika_XrayEnergy="+num2str(XrayEnergy)+";"
		OldNote+= "Nika_PixelSizeX="+num2str(PixelSizeX)+";"
		OldNote+= "Nika_PixelSizeY="+num2str(PixelSizeY)+";"
		OldNote+= "Nika_HorizontalTilt="+num2str(HorizontalTilt)+";"
		OldNote+= "Nika_VerticalTilt="+num2str(VerticalTilt)+";"
		OldNote+= "Nika_BeamCenterX="+num2str(BeamCenterX)+";"
		OldNote+= "Nika_BeamCenterY="+num2str(BeamCenterY)+";"
		OldNote+= "Nika_BeamSizeX="+num2str(BeamSizeX)+";"
		OldNote+= "Nika_BeamSizeY="+num2str(BeamSizeY)+";"
		OldNote+= "CalibrationFormula="+CalibrationFormula+";"
		if(UseSampleThickness)
			OldNote+= "SampleThickness="+num2str(SampleThickness)+";"
		endif
		if(UseSampleTransmission)
			OldNote+= "SampleTransmission="+num2str(SampleTransmission)+";"
		endif
		if(UseCorrectionFactor)
			OldNote+= "CorrectionFactor="+num2str(CorrectionFactor)+";"
		endif
		if(UseSubtractFixedOffset)
			OldNote+= "SubtractFixedOffset="+num2str(SubtractFixedOffset)+";"
		endif
		if(UseSampleMeasTime)
			OldNote+= "SampleMeasurementTime="+num2str(SampleMeasurementTime)+";"
		endif
		if(UseEmptyMeasTime)
			OldNote+= "EmptyMeasurementTime="+num2str(EmptyMeasurementTime)+";"
		endif
		if(UseI0ToCalibrate)
			OldNote+= "SampleI0="+num2str(SampleI0)+";"
			OldNote+= "EmptyI0="+num2str(EmptyI0)+";"
		endif
		if(UseDarkMeasTime)
			OldNote+= "BackgroundMeasTime="+num2str(BackgroundMeasTime)+";"
		endif
		if(UsePixelSensitivity)
			OldNote+= "UsedPixelsSensitivity="+num2str(UsePixelSensitivity)+";"
		endif
		if(UseMonitorForEF)
			OldNote+= "UseMonitorForEF="+num2str(UseMonitorForEF)+";"
		endif
	
		SVAR CurrentDarkFieldName=root:Packages:Convert2Dto1D:CurrentDarkFieldName
		SVAR CurrentEmptyName=root:Packages:Convert2Dto1D:CurrentEmptyName	
		if(UseDarkField)
			OldNote+= "CurrentDarkFieldName="+(CurrentDarkFieldName)+";"
		endif
		if(UseEmptyField)
			OldNote+= "CurrentEmptyName="+(CurrentEmptyName)+";"
		endif
	else
		OldNote+= "CalibrationFormula="+"1"+";"
	endif
	SVAR CurrentMaskFileName=root:Packages:Convert2Dto1D:CurrentMaskFileName
	if(UseMask)
		OldNote+= "CurrentMaskFileName="+(CurrentMaskFileName)+";"
	endif
	NVAR UseSolidAngle=root:Packages:Convert2Dto1D:UseSolidAngle
	if(UseSolidAngle)
		OldNote+= "SolidAngleCorrection=Done"+";"
	endif
	
	note /K Calibrated2DDataSet
	note Calibrated2DDataSet, OldNote
	KillWaves/Z tempEmptyField, tempDataWv
	setDataFolder OldDf
end

//*******************************************************************************************************************************************
//*******************************************************************************************************************************************

Function/WAVE SolidAngleCorrection(DataWave)
	Wave DataWave
	Wave/Z Param_Values = root:Param_Values
	
	if(WaveExists(Param_Values) == 0)
		print "Making Needed Paramater Waves with Default Values"
		setDataFolder root:
		Make/T=2 Param_Names = {"Distance from axis of rotation to sample (mm)", "Sample - Detector distance (mm)"}
		Make/N=2 Param_Values = {26,23}

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
	
	variable ROI_Width = str2num(Stringbykey("Camera ROI Width",header,"=",";"))
	variable ROI_Height = str2num(Stringbykey("Camera ROI Height",header,"=",";"))
	variable ROI_X_Bin = str2num(Stringbykey("Camera ROI X Bin",header,"=",";"))
	variable ROI_Y_Bin = str2num(Stringbykey("Camera ROI Y Bin",header,"=",";"))
	
	variable num_pixels_X = ROI_Width / ROI_X_Bin
	variable num_pixels_Y = ROI_Height / ROI_Y_Bin
	
	NVAR PixelSizeX = root:Packages:Convert2Dto1D:PixelSizeX							//in millimeters
	NVAR PixelSizeY = root:Packages:Convert2Dto1D:PixelSizeY							//in millimeters
	NVAR beamCenterX=root:Packages:Convert2Dto1D:beamCenterX
	NVAR beamCenterY=root:Packages:Convert2Dto1D:beamCenterY

	Variable C = BeamCenterY*PixelSizeY	//Distance from the beam center in normal incidence to bottom of the CCD (mm)

	Duplicate/O DataWave, SolidAngleCorrectionMap
	
	Make/Free/N=(num_pixels_Y) a
	Make/Free/N=(num_pixels_X) b
	
	a = (p*PixelSizeY) + (d1 + d2 - (d2 + d1*(1-cos(phi))))*tan(phi) - C //Place Holder variable for y-axis 1D WAVE
	b = (p - beamCenterX ) * PixelSizeX //Place Holder variable for x-axis 1D WAVE

	//Turn Everything into Matrices so elementwise low level matrix operations can be done.
	MatrixOP/FREE y1 = rowRepeat(a,(num_pixels_Y))
	MatrixOP/FREE y2 = y1 + PixelSizeY
	MatrixOP/FREE x1 = colRepeat(b,(num_pixels_X))
	MatrixOP/FREE x2 = x1 + PixelSizeX
	MatrixOP/FREE z1 = const((num_pixels_Y),(num_pixels_X),(d2 + d1*(1-cos(phi))))
			// z is the distance from the source perpendicular to the Cartesian plane the CCD is in.
	
	MatrixOP/FREE omega = atan(x2*y2/(z1*sqrt(sq(x2)+sq(y2)+sq(z1)))) - atan(x2*y1/(z1*sqrt(sq(x2)+sq(y1)+sq(z1)))) - atan(x1*y2/(z1*sqrt(sq(x1)+sq(y2)+sq(z1))))	+ atan(x1*y1/(z1*sqrt(sq(x1)+sq(y1)+sq(z1))))
	// omega is the solid angle of a pixel
	MatrixOP/FREE omega_0 = 4*atan(PixelSizeX*PixelSizeY/(2*z1*sqrt(sq(PixelSizeX)+sq(PixelSizeY)+4*sq(z1))))

	MatrixOP/O SolidAngleCorrectionMap = omega_0 / omega //This is a normalization of the solid angle to virtual pixel residing at the termination of 'z'
	
	return SolidAngleCorrectionMap

end