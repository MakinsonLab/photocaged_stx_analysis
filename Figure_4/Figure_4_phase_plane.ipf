//Uses abf import function from NeuroMatic Toolkit written by Rothman JS and Silver RA
//NeuroMatic: An Integrated Open-Source Software Toolkit for Acquisition, 
//Analysis and Simulation of Electrophysiological Data. 
//Front Neuroinform. 2018 Apr 4;12:14. doi: 10.3389/fninf.2018.00014.


Function ListFilesSTX()
SVAR  gCustomPath
Variable NumFiles, i
String FileList, abf_files, Current_file,Current_abf,Parent_file,Cell_file

NewDataFolder/S/O root:AP_numbers

Make/O/N=0 NumberAPsTbl
Make/O/T/N=0 FileNameTbl
Make/O/N=0 SweepNumberTbl
Make/O/T/N=0 DateTbl
Make/O/N=0 TimeMSTbl
Display/N=PP_stx
Display/N=AP_stx
Display/N=D_stx
Display/N=SD_stx


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


AP_Calc(Current_file,Current_abf)

endfor

yay_rainbow("PP_stx")
yay_rainbow("AP_stx")
yay_rainbow("D_stx")
yay_rainbow("SD_stx")

Edit FileNameTbl,SweepNumberTbl,NumberAPsTbl,TimeMSTbl,DateTbl

ModifyTable/W=Table0 format(TimeMSTbl)=1
ModifyTable/W=Table0 format(DateTbl)=1

string save_file_path = gCustomPath+Parent_file+"_"+Cell_file+"_stx.txt"

//print save_file_path

SaveTableCopy/O/W=Table0/T=1 as save_file_path
KillWindow Table0


NewLayout/K=1/N=CagedSweeps as "Phase Plane Plot for caged Stx"
AppendLayoutObject/W=CagedSweeps graph PP_stx
AppendLayoutObject/W=CagedSweeps graph AP_stx
AppendLayoutObject/W=CagedSweeps graph D_stx
AppendLayoutObject/W=CagedSweeps graph SD_stx

Execute("Tile/O=0")


string pdf_file_path = gCustomPath+Parent_file+"_"+Cell_file+"_stx.pdf"

SavePict/E=-8/O/WIN=CagedSweeps as pdf_file_path




End


Function AP_Calc(Current_file,Current_abf)
String Current_file,Current_abf

DFREF CurrDF
variable j,k, Number_APs,SweepNumber,FileMod,NumberAPsInSweep
string TrcPotNm, TrcCmdCurrNm, Start_date

setDataFolder root:AP_numbers
Wave NumberAPsTbl, SweepNumberTbl,TimeMSTbl//,DateTbl
Wave/T FileNameTbl,DateTbl
setDataFolder root:





NMImportFile( "new" ,Current_file)
DoWindow /K ImportPanel
SVAR FileName
String FileNameNoExt = RemoveEnding(FileName , ".abf")

setdataFolder ABFHeader
NVAR uFileStartTimeMS
NVAR uFileStartDate
//print uFileStartDate

Start_date = num2istr(uFileStartDate)

//print Start_date

setdataFolder ::
 

NVAR FileDateTime
NVAR Num_traces = NumWaves



variable stepstarttime = 45
variable stependtime = 1045
string AP_name
variable beendone = 0

j = 0


//Each sweep
for(j=0;j<Num_traces;j+=1)
	NumberAPsInSweep = 0

	TrcPotNm = "RecordA" +num2str(j)
 	Wave TracePot = $TrcPotNm
 	TrcCmdCurrNm = "RecordB" +num2str(j)
 	Wave TraceCmdCurr = $TrcCmdCurrNm
 
 	String  TrcPotDiffNm = TrcPotNm+"_dif"
 	string AP_times_Nm = "RecordA" +num2str(j)+"_AP_times"

 Smooth 20, TracePot
 	Differentiate TracePot/D=$TrcPotDiffNm
 	

 //Find action potentials within sweep
 	FindLevels/EDGE=1/Q/R=(stepstarttime,stependtime)/D=$AP_times_Nm $TrcPotDiffNm, 5
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
			APEndA = (AP_times[k])+10	
			
			print ("k = "+num2str(k))
			print ("j = "+num2str(j))			
		
			
			if (j==40)
			print("only me")
				
				
				AP_name = FileName+"_S_"+num2str(k)
				Duplicate/O/R=(APStartA,APEndA) TracePot $AP_name
				String AP_name_diff = FileNameNoExt+"_S_"+num2str(k)+"_dif"
				String AP_name_sec_diff = FileNameNoExt+"_S_"+num2str(k)+"_dif_dif"				
				Differentiate $AP_name/D=$AP_name_diff
				Differentiate $AP_name_diff/D=$AP_name_sec_diff
				
				
				
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

			
		Endfor
		k=0

endif


	

 SweepNumber = j+1
 
 CurrDF = GetDataFolderDFR()
 
 setDataFolder root:AP_numbers
 
// print num2str(SweepNumber)
// print Current_abf
 	
 	InsertPoints inf,1, NumberAPsTbl
	NumberAPsTbl[inf] =  NumberAPsInSweep

	InsertPoints inf,1, SweepNumberTbl
	SweepNumberTbl[inf] =  SweepNumber
	
	InsertPoints inf,1, FileNameTbl
	FileNameTbl[inf] =  Current_abf
	
	InsertPoints inf,1, TimeMSTbl
	TimeMSTbl[inf] =  uFileStartTimeMS
	
	InsertPoints inf,1, DateTbl
	DateTbl[inf] =  Start_date
	
	
	
	

 
 setDataFolder CurrDF
 
	
 
endfor
 
 
 

End


Function yay_rainbow(thebigname)
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


