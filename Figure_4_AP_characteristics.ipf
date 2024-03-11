 //Uses abf import function from NeuroMatic Toolkit written by Rothman JS and Silver RA
//NeuroMatic: An Integrated Open-Source Software Toolkit for Acquisition, 
//Analysis and Simulation of Electrophysiological Data. 
//Front Neuroinform. 2018 Apr 4;12:14. doi: 10.3389/fninf.2018.00014.

Function Advanced_analysis(current_abf_file,RecInfo)
String current_abf_file
Variable RecInfo

setdatafolder root:
Wave Sweep
Wave/T Trace_name,Condition,Slice

Wave Threshold_AP_amp_Tbl,AP_amp_Tbl,dif_amp_Tbl,p1_p2_Tbl,p1_amp_Tbl,Sweep_number_Tbl,p2_latency_Tbl
Wave/T Condition_Tbl,Trace_name_1_Tbl,Trace_name_2_Tbl,Slice_Tbl



Variable SweepNumber = Sweep[RecInfo]
String Condition_nm = Condition[RecInfo]
String Trace_name_nm = Trace_name[RecInfo]
String Slice_nm = Slice[RecInfo]
print(num2str(RecInfo))
print(current_abf_file)
//Correct file path 
current_abf_file = ReplaceString("::", current_abf_file, ":")
print(current_abf_file)

DFREF CurrDF
variable j,k, Number_APs,FileMod,NumberAPsInSweep
string TrcPotNm, TrcCmdCurrNm, Start_date

//setDataFolder root:AP_numbers
//Wave NumberAPsTbl, SweepNumberTbl,TimeMSTbl//,DateTbl
//Wave/T FileNameTbl,DateTbl
setDataFolder root:


NMImportFile( "new" ,current_abf_file)
DoWindow /K ImportPanel
SVAR FileName
String FileNameNoExt = RemoveEnding(FileName , ".abf")

NVAR Num_traces = NumWaves



variable stepstarttime = 110
variable stependtime = 615
string AP_name
variable beendone = 0



//Select sweep 

TrcPotNm = "RecordA" +num2str(SweepNumber-1)
Wave TracePot = $TrcPotNm

Display $TrcPotNm

String  TrcPotDiffNm = TrcPotNm+"_dif"
//String  TrcPotSecDiffNm = TrcPotNm+"_dif_dif"

String AP_times_Nm = "RecordA" +num2str(j)+"_AP_times"

Smooth 20, TracePot
Differentiate TracePot/D=$TrcPotDiffNm
 	
//Find action potentials within sweep
 	FindLevels/EDGE=1/Q/R=(stepstarttime,stependtime)/D=$AP_times_Nm $TrcPotDiffNm, 10
 	Print("***********************************")
 	
 	PRint("File ID: "+FileNameNoExt)
 	PRint("Sweep Number: "+TrcPotNm)
 	PRint("Number of APs: "+num2str(V_LevelsFound))
	if(V_LevelsFound>0)
	
		Wave AP_times = $AP_times_Nm
		PRint("AP Number Wave: "+AP_times_Nm)
		print("Number APs detected: "+num2str(numpnts(AP_times)))
		
		Variable APStartA, APEndA, AP_VmaxA
		
		//First AP in sweep analysis		

			APStartA =(AP_times[0])-4
			APEndA = (AP_times[0])+10			
				
			//Create snips of 1st AP in selected sweep
			AP_name = FileNameNoExt+"_S_"+num2str(SweepNumber)
			Duplicate/O/R=(APStartA,APEndA) TracePot $AP_name
			Wave AP_snip = $AP_name
			
			//Differentiate AP
			String AP_name_diff = FileNameNoExt+"_S"+num2str(SweepNumber)+"_dif"
			String AP_name_sec_diff = FileNameNoExt+"_S_"+num2str(SweepNumber)+"_dif_dif"				
			Differentiate $AP_name/D=$AP_name_diff
			Differentiate $AP_name_diff/D=$AP_name_sec_diff
			
			//Measurements
			//Threshold
			
			Display AP_snip
			FindLevel/EDGE=1/P/Q $AP_name_diff, 5
			Variable Vthresholdpoint = V_LevelX
			Variable V_threshold = AP_snip[round(V_LevelX)]
			//print(num2str(Vthresholdpoint))
			print("Voltage threshold: "+num2str(V_threshold))
			
			//Amplitude
			Wavestats/Q $AP_name
			Variable Amplitude = V_max - V_threshold
			print("AP amplitude: "+num2str(Amplitude))
			
			//Max
			Wavestats/Q $AP_name_diff
			Variable Maxdif = V_max
			print("Diff amplitude: "+num2str(Maxdif))
			
			string man_plot_name = "Plot_"+FileNameNoExt+"_S_"+num2str(SweepNumber)
			Wave DoubleDif = $AP_name_sec_diff
			Display/N=$man_plot_name $AP_name_sec_diff
			DoWindow $man_plot_name			
			ShowInfo
	

		if (UserCursorAdjust(man_plot_name,30) != 0)
			return -1
		endif

		if (strlen(CsrWave(A))>0 && strlen(CsrWave(B))>0)	// Cursors are on trace?
				
				
				variable hump_one = pcsr(A)
				variable hump_two = pcsr(B)
				
				print "-------------------------------------------"
				
				variable double_diff_peak_one = DoubleDif[hump_one]
				variable double_diff_peak_two = DoubleDif[hump_two]
				variable peak_ratio = double_diff_peak_one/double_diff_peak_two
				
				print("Peak one amplitude: "+num2str(DoubleDif[hump_one]))
				print("Peak two amplitude: "+num2str(DoubleDif[hump_two]))
				
				variable double_diff_peak_time_one = IndexToScale(DoubleDif,hump_one,0)
				variable double_diff_peak_time_two = IndexToScale(DoubleDif,hump_two,0)
				
				variable double_diff_latency = double_diff_peak_time_two - double_diff_peak_time_one
				print("Peak one time: "+num2str(double_diff_peak_time_one))
				print("Peak two time: "+num2str(double_diff_peak_time_two))
				
				
				print("lATENCY: "+NUM2STR(double_diff_latency))
				
				
				
				
				
				
		endif
		
	
		
	InsertPoints inf,1, Threshold_AP_amp_Tbl
	Threshold_AP_amp_Tbl[inf] =  V_threshold

	InsertPoints inf,1, AP_amp_Tbl
	AP_amp_Tbl[inf] =  Amplitude
	
	InsertPoints inf,1, dif_amp_Tbl
	dif_amp_Tbl[inf] = Maxdif
	 
	//InsertPoints inf,1, dif_dif_Tbl
	//dif_dif_Tbl[inf] = double_diff_peak_one 
	
	InsertPoints inf,1, p1_p2_Tbl
	p1_p2_Tbl[inf] =  peak_ratio
	
	InsertPoints inf,1, p1_amp_Tbl
	p1_amp_Tbl[inf] =  double_diff_peak_one 
	
	InsertPoints inf,1, Condition_Tbl
	Condition_Tbl[inf] =  Condition_nm
	
	InsertPoints inf,1, Sweep_number_Tbl
	Sweep_number_Tbl[inf] =  SweepNumber
	
	InsertPoints inf,1, Trace_name_1_Tbl
	Trace_name_1_Tbl[inf] = Trace_name_nm 
	
	InsertPoints inf,1, Trace_name_2_Tbl
	Trace_name_2_Tbl[inf] =  current_abf_file
	
	InsertPoints inf,1, Slice_Tbl
	Slice_Tbl[inf] = Slice_nm
	
	InsertPoints inf,1, p2_latency_Tbl
	p2_latency_Tbl[inf] = double_diff_latency
	
			
			
			
			
			
			
		endif
			
End				
				
				
				AppendToGraph/W=PP_stx $AP_name_diff vs $AP_name
				
				//
				String  Fig_AP = FileName+"_fig_S_"+num2str(k)
				Variable Fig_AP_start =(AP_times[k])-0.25
				Variable Fig_AP_end = (AP_times[k])+5	
				Duplicate/O/R=(Fig_AP_start,Fig_AP_end) TracePot $Fig_AP
				Setscale/P x,0,deltax($Fig_AP), $Fig_AP	
				String AP_fig_diff = FileNameNoExt+"_S_"+num2str(k)+"_fig_dif"
				String AP_fig_sec_diff = FileNameNoExt+"_S_"+num2str(k)+"_fig_dif_dif"
				Differentiate $Fig_AP/D=$AP_fig_diff
				Differentiate $AP_fig_diff/D=$AP_fig_sec_diff		
						
				AppendToGraph/W=AP_stx $Fig_AP
				AppendToGraph/W=D_stx $AP_fig_diff
				AppendToGraph/W=SD_stx $AP_fig_sec_diff	
				
				
				
				
						
			
			endif

			Wavestats/Q/R=(APStartA,APEndA)Tracepot

			AP_VmaxA = V_max

			//if(V_max>0)
			//	NumberAPsInSweep = NumberAPsInSweep+1				
			//endif

		Endfor
		k=0

endif


	

 SweepNumber = j+1
 
 CurrDF = GetDataFolderDFR()
 
// setDataFolder root:AP_numbers
 
// print num2str(SweepNumber)
// print Current_abf
 	
 	



	
	
	
setDataFolder CurrDF
 
	
 
endfor
 
 
 

End


#pragma rtGlobals = 1

Function yewq(thebigname)
string thebigname
string PlotNames  
variable numTraces
ColorTab2Wave BlueRedGreen

//TraceNameList(graphNameStr, separatorStr, optionsFlag)

PlotNames = (TraceNameList(thebigname, ";", 1))

print(PlotNames)
numTraces = itemsinList(PlotNames)
print (num2str(numTraces))


	if (numTraces <= 0)
		return -1
	endif
	
   Variable denominator= numTraces-1
   if( denominator < 1 )
       denominator= 1    // avoid divide by zero, use just the first color for 1 trace
   endif

Wave rgb = M_colors
	Variable numRows= DimSize(rgb,0)
	Variable red, green, blue
	Variable i, index
	for(i=0; i<numTraces; i+=1)
		index = round(i/denominator * (numRows-1))	// spread entire color range over all traces.
		ModifyGraph/W=$thebigname rgb[i]=(rgb[index][0], rgb[index][1], rgb[index][2])
	endfor



End

//-----------------------------------------------------------

Function UserCursorAdjust(graphName,autoAbortSecs)
	String graphName
	Variable autoAbortSecs

	DoWindow/F $graphName							// Bring graph to front
	if (V_Flag == 0)									// Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif

	NewPanel /K=2 /W=(187,368,437,531) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor					// Set to an unlikely name
	AutoPositionWindow/E/M=1/R=$graphName			// Put panel near the graph

	DrawText 21,20,"Adjust the cursors and then"
	DrawText 21,40,"Click Continue."
	Button button0,pos={80,58},size={92,20},title="Continue"
	Button button0,proc=UserCursorAdjust_ContButtonProc
	Variable didAbort= 0
	if( autoAbortSecs == 0 )
		PauseForUser tmp_PauseforCursor,$graphName
	else
		SetDrawEnv textyjust= 1
		DrawText 162,103,"sec"
		SetVariable sv0,pos={48,97},size={107,15},title="Aborting in "
		SetVariable sv0,limits={-inf,inf,0},value= _NUM:10
		Variable td= 10,newTd
		Variable t0= ticks
		Do
			newTd= autoAbortSecs - round((ticks-t0)/60)
			if( td != newTd )
				td= newTd
				SetVariable sv0,value= _NUM:newTd,win=tmp_PauseforCursor
				if( td <= 10 )
					SetVariable sv0,valueColor= (65535,0,0),win=tmp_PauseforCursor
				endif
			endif
			if( td <= 0 )
				DoWindow/K tmp_PauseforCursor
				didAbort= 1
				break
			endif
				
			PauseForUser/C tmp_PauseforCursor,$graphName
		while(V_flag)
	endif
	return didAbort
End

Function UserCursorAdjust_ContButtonProc(ctrlName) : ButtonControl
	String ctrlName

	DoWindow/K tmp_PauseforCursor				// Kill self
End

Function Demo()
	
	Wave Record0
	
	
	
		Display/N=yup RecordA0
		ShowInfo
	

	if (UserCursorAdjust("yup",10) != 0)
		return -1
	endif

	if (strlen(CsrWave(A))>0 && strlen(CsrWave(B))>0)	// Cursors are on trace?
		print(pcsr(A))
		print(pcsr(B))
		Print IndexToScale(w4D,1,0)
	endif
End
