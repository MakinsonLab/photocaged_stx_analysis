 //Uses abf import function from NeuroMatic Toolkit written by Rothman JS and Silver RA
//NeuroMatic: An Integrated Open-Source Software Toolkit for Acquisition, 
//Analysis and Simulation of Electrophysiological Data. 
//Front Neuroinform. 2018 Apr 4;12:14. doi: 10.3389/fninf.2018.00014.




#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Function gap_free_AP_analysis()
string  gCustomPath = "C:/Users/"
Variable NumFiles, i
String FileList, abf_files, Current_file,Current_abf,Parent_file,Cell_file


LoadWave/A/W/J/N/O/K=2

Wave/T stx_trace_name
Wave/T UV_on



///NewDataFolder/S/O root:AP_numbers

SetDataFolder root:

Make/O/N=0 AP_number_table
Make/O/T/N=0 FileNameTbl
Make/O/N=0 Amplitude_table
Make/O/N=0 Threshold_table
Make/O/N=0 AP_time_table
Make/O/N=0 UV_on_time_table
Make/O/N=0 Time_after_UV_table
Make/O/N=0 Peak_dv_dt


NewPath/M="Cell Folder"/O DirectoryPath

PathInfo DirectoryPath

Parent_file = ParseFilePath(0, S_path, ":", 1, 1)
//Parent_file = ParseFilePath(0, S_path, ":", 1, 2)

Cell_file = ParseFilePath(0, S_path, ":", 1, 0)

FileList = IndexedFile(DirectoryPath, -1, "????")

abf_files =GrepList(FileList,".*\.abf" )
abf_files = SortList(abf_files)

print(abf_files)

NumFiles = ItemsInList (abf_files)




for (i=0;i<NumFiles;i+=1)

Current_abf = StringFromList(i,abf_files)

Current_file = S_path+Current_abf

String match_name = RemoveEnding(Current_abf , ".abf")



FindValue/TEXT=match_name/TXOP=2 stx_trace_name


string uv_time_string = UV_on[V_row]
string stx_string = stx_trace_name[V_row]

print("---------------------------------------------------")
print(match_name)
print(uv_time_string)
print(stx_string)
print("---------------------------------------------------")


variable UV_time = str2num(uv_time_string) 



AP_GF(Current_file,Current_abf,UV_time)

endfor

Edit AP_number_table,FileNameTbl,Amplitude_table,Threshold_table,AP_time_table,UV_on_time_table,Time_after_UV_table,Peak_dv_dt

string save_file_path = gCustomPath+Parent_file+"_"+Cell_file+"_gap_stx.txt"

print save_file_path

SaveTableCopy/O/W=Table0/T=1 as save_file_path
KillWindow Table0

End


Function AP_GF(Current_file,Current_abf,UV_time)



String Current_file,Current_abf
variable UV_time

DFREF CurrDF
variable j,k, Number_APs,SweepNumber,FileMod,NumberAPsInSweep
string  gCustomPath = "C:/Users/"

string TrcPotNm, TrcCmdCurrNm

setDataFolder root:
Wave/T  FileNameTbl
Wave AP_number_table,Amplitude_table,Threshold_table,AP_time_table,UV_on_time_table,Peak_dv_dt


string promptstring = "Light Info"+ Current_abf



//Prompt UV_time "UV on time"
//DoPrompt  promptstring ,UV_time

NMImportFile( "new" ,Current_file)
DoWindow /K ImportPanel
SVAR FileName
String FileNameNoExt = RemoveEnding(FileName , ".abf")

String Plot_nm = "G_"+FileNameNoExt+"_plot"
Display/N=$Plot_nm


 TrcPotNm = "RecordA" +num2str(0)
 Wave TracePot = $TrcPotNm
 	
 
String  TrcPotDiffNm = TrcPotNm+"_dif"
string AP_times_Nm = "RecordA" +num2str(j)+"_AP_times"


	Smooth 5, TracePot 
  	Differentiate TracePot/D=$TrcPotDiffNm
 	

 
 	FindLevels/EDGE=1/Q/D=$AP_times_Nm $TrcPotDiffNm, 30
 	Print("***********************************")
 	
 	PRint("File ID: "+FileNameNoExt)
 	PRint("Sweep Number: "+TrcPotNm)
 	PRint("Number of APs: "+num2str(V_LevelsFound))
	if(V_LevelsFound>0)
	
		Wave AP_times = $AP_times_Nm
		PRint("AP Number Wave: "+AP_times_Nm)
		print("Number APs detected: "+num2str(numpnts(AP_times)))		
		Variable APStartA, APEndA, AP_VmaxA
		For(k=0;k<numpnts(AP_times);k+=1)	

			APStartA =(AP_times[k])-4
			APEndA = (AP_times[k])+6	
			
			print ("k = "+num2str(k))
			
			//Create snips of 1st AP in selected sweep
			String AP_name = FileNameNoExt+"_"+num2str(k)
			Duplicate/O/R=(APStartA,APEndA) TracePot $AP_name
			Wave AP_snip = $AP_name
			
			//Differentiate AP
			String AP_name_diff = FileNameNoExt+"_S"+num2str(SweepNumber)+"_dif"					
			Differentiate $AP_name/D=$AP_name_diff
			Wavestats/Q $AP_name_diff
			variable dvdt = V_max
			print("max dv/dt: "+num2str(dvdt)	)			
			
			//Threshold			
			FindLevel/EDGE=1/P/Q $AP_name_diff, 5
			Variable Vthresholdpoint = V_LevelX
			Variable V_threshold = AP_snip[round(V_LevelX)]			
			print("Voltage threshold: "+num2str(V_threshold))
			
			//Amplitude
			Wavestats/Q $AP_name
			Variable Amplitude = V_max - V_threshold
			print("AP amplitude: "+num2str(Amplitude))	
			
			AppendToGraph/W=$Plot_nm AP_snip
			
			
		DFREF CurrDF = getDataFolderDFR()
		setdataFolder root:

		
		InsertPoints inf,1, FileNameTbl
		FileNameTbl[inf] =  Current_abf
		
		InsertPoints inf,1, AP_number_table
		AP_number_table[inf] =  k+1
	
		InsertPoints inf,1, Amplitude_table
		Amplitude_table[inf] =  Amplitude
	
		InsertPoints inf,1, Threshold_table
		Threshold_table[inf] =  V_threshold
		
		InsertPoints inf,1, AP_time_table
		AP_time_table[inf] =  AP_times[k]
		
		InsertPoints inf,1, UV_on_time_table
		UV_on_time_table[inf] =  UV_time*1000
		
		InsertPoints inf,1, Peak_dv_dt
		Peak_dv_dt[inf] =  dvdt
		
		
		setDataFolder CurrDF	
			

		Endfor
		
	variable plot_start,plot_stop
	plot_start = (UV_time-5)*1000
	plot_stop = (UV_time+10)*1000
	
	duplicate/O TracePot TracePlot_color_index, red, green, blue
	
	
	
	variable point_of_UV = x2pnt(TracePot, UV_time*1000 )
	
	variable idx
	
	for(idx=0;idx<numpnts(TracePot);idx+=1)
	
	red[idx]	= 65535
	green[idx] = 0
	blue[idx] = 0
	
	if(idx > point_of_UV)
	
	red[idx]	= 0	
	blue[idx] = 65535
	
	endif	
	
	endfor
	

	
	Make/O/N=(numpnts(red),3) color_wave
	color_wave[][0] = red[p] // Store into all rows, column 0
	color_wave[][1] = green[p] // Store into all rows, column 1
	color_wave[][2] = blue[p]	
	
	
	Make/O/N=(numpnts(TracePot)) fzwave = p
	Display/N=colored_plot  TracePot
	ModifyGraph/W=colored_plot zColor(RecordA0)={color_wave,*,*,directRGB}
	SetAxis/W=colored_plot bottom, plot_start,plot_stop
	
	
	//string save_file_path = gCustomPath+cell_name_no_ext+"_gap_stx.txt"

	//SaveTableCopy/O/W=Table0/T=1 as save_file_path
	//KillWindow Table0

	string save_trace_path = gCustomPath+FileNameNoExt+"_gap_stx.pdf"
	SavePict/E=-8/O/WIN=colored_plot as save_trace_path
	
	KillWindow colored_plot

		
		
	
	endif
	
	
	
	

 
 
End