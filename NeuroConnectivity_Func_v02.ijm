/*	NeuroConnectivity_Func_v2.ijm
	*************4**********
	Author: 			Winnok H. De Vos
	Modified by: 		Marlies Verschuuren
	Date Created: 		December 17th, 2010
	Date Last Modified:	Januari 23th, 2024 
 	
 	Description: 
 	Macro Set dedicated to the analysis of calcium fluxes in whole fields or cellular ROIs within the field of view
 	
 	Change Log
 	* Solve bug single image analysis if no additional marker is used (250108-Marlies Verschuuren-v2)
 	_________________________________________________________________
*/	

/*
 	***********************

	Variable initiation

	***********************
*/

//	String variables
var cell_threshold					= "Fixed";								//	threshold method for segmentation of cells 
var cell_method						= "Stardist";								// 	segmentation method
var dir								= "";										//	directory
var file_list						= "";
var log_path						= "";										//	path for the log file
var marker_identifier				= "_MarkerRoi"								// 	string to identify file with additional marker
var micron							= getInfo("micrometer.abbreviation");		// 	micro symbol
var output_dir						= "";										//	output directory
var roi_type						= "Additional Marker"									//	ROI(s) to be measured
var suffix_time						= ".tif";									//	suffix for specifying the time file type
var suffix_marker					= ".tif";									//	suffix for specifying the file type

//	Numeric variables
var cell_area_min					= 50;										//  approximate min area of cells 
var cell_area_max					= 700;										//  approximate max area of cells 
var cell_filter_scale				= 0;										//	radius for cell smoothing
var cell_fixed_threshold_value		= 1000;										//	fixed maximum threshold for cell segmentation (if auto doesn't work well);
var image_height					= 1000;										//	image height
var image_width						= 1000;										//	image width
var marker_scale					= 0.5;										// 	scale factor marker image
var pixel_size						= 1.434;									//	pixel size (µm)
var prominence						= 25										//	prominence for max finding and particle dectetion
var stardist_probability			= 0.3
var stardist_overlap				= 0.15

//	Array variables
var file_types 						= newArray(".tif",".tiff",".nd2",".ids",".jpg",".mvd2",".czi");		
var images 							= newArray(0);
var list							= newArray(0);
var roi_types						= newArray("Full FOV","Calcium Signal","Additional Marker");
var prefixes 						= newArray(0);
var segmentation_methods			= newArray("Thresholding","Stardist");
var threshold_methods				= getList("threshold.methods");	
var threshold_methods				= Array.concat(threshold_methods,"Fixed");	

/*
 	***********************

		Macros

	***********************
*/

macro "Setup Action Tool - C888 T5f16S"
{
	setup();
}

macro "Analyse Single Image Action Tool - C888 T5f161"
{
	erase(0);
	setBatchMode(true);
	id = getImageID;
	calibrateImage(id);
	dir=getDirectory("image");
	title=getTitle();
	prefix=substring(title,0,lastIndexOf(title,suffix_time));
	measureStack(id,roi_type, prefix);
	selectImage(id);
	getDimensions(width, height, channels, slices, frames);
	if(roi_type=="Additional Marker"){
		path = dir+prefix+marker_identifier+suffix_marker;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		idm=getImageID();
		run("Scale...", "x="+marker_scale+" y="+marker_scale+" width="+width+" height="+height+" interpolation=Bilinear average create");
		idmDup=getImageID();
		decalibrateImage(idmDup);
		calibrateImage(idmDup);
		roiManager("show all without labels");
		selectImage(idm); close();
	}
	
	setBatchMode("exit and display");
	run("Tile");
}

macro "Batch Analysis Action Tool - C888 T5f16#"
{	
	erase(0);
	setOptions();
	setDirectory();
	prefixes = scanFiles();
	image_nr = prefixes.length;
	start = getTime();
	setBatchMode(true);
	for(i = 0; i < image_nr; i++)
	{
		prefix = prefixes[i];
		if(roi_type=="Additional Marker" && prefix.indexOf(marker_identifier)>0){
		}else{
			file = prefix+suffix_time;
			print(i+1,"/",image_nr,":",prefix);
			path = dir+file;
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
			id 	= getImageID;
			calibrateImage(id);
			measureStack(id,roi_type, prefix);
			selectWindow("Results");
			saveAs("Measurements", output_dir+prefix+"_results.txt");
			roiManager("Save",output_dir+prefix+"_roi_set.zip"); 
			selectImage(id); 
			close();
			run("Clear Results");
			roiManager("reset");
			run("Collect Garbage");
		}
	}
	print((getTime()-start)/1000,"sec");
	if(isOpen("Log"))
	{
		selectWindow("Log");
		saveAs("txt",log_path);
	}
	print("Analysis Done");
	setBatchMode("exit and display");
}

macro "Toggle Overlay Action Tool - Caaa O11ee"
{
	toggleOverlay();
}

macro "[t] Toggle Overlay"
{
	toggleOverlay();
}

macro "Verification Stack Action Tool - C888 T5f16V"
{
	erase(1);
	setBatchMode(true);
	setDirectory();
	prefixes = scanFiles();
	names=newArray(0);
	if(roi_type=="Additional Marker")
	{
		for(i=0;i<prefixes.length;i++)
		{
			if(prefixes[i].indexOf(marker_identifier)<0)
			{
				names=Array.concat(names,prefixes[i]);
			}
		}
	}else{
		names = prefixes;
	}
	createOverlay(names);
	setBatchMode("exit and display");
	run("Channels Tool... ");
	run("Tile");
	run("Synchronize Windows");
}

/*
 	***********************

		Functions

	***********************
*/

function setOptions()
{
	run("Options...", "iterations=1 count=1");
	run("Colors...", "foreground=white nuclei_background=black selection=yellow");
	run("Overlay Options...", "stroke=red width=1 fill=none");
	setBackgroundColor(0, 0, 0);
	setForegroundColor(255,255,255);
	run("Set Measurements...", "mean redirect=None decimal=3");
}

function getMoment()
{
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
     return TimeString;
}

function erase(all)
{
	if(all){run("Close All");}
	print("\\Clear");
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
}

function setDirectory()
{
	dir = getDirectory("Choose a Source Directory");
	file_list = getFileList(dir);
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir))File.makeDirectory(output_dir);
	log_path = output_dir+"Log.txt";
}

function scanFiles()
{
	prefixes = newArray(0);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,suffix_time) && indexOf(path,"flatfield")<0)
		{
			print(path);
			prefixes = Array.concat(prefixes,substring(file_list[i],0,lastIndexOf(file_list[i],suffix_time)));			
		}
	}
	return prefixes;
}

function setup()
{
	setOptions();
	Dialog.createNonBlocking("Calcium Analysis Settings");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("General:\n", 13, "#f26852");	
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Time Serie: Image Type", file_types, suffix_time);
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Pixel Size", pixel_size, 3, 5, micron+"");
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("ROI Type", roi_types, roi_type);
	Dialog.setInsets(0,0,0);
	Dialog.addString("> Marker: identifier", marker_identifier, 15);
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("> Marker: Image Type", file_types, suffix_marker);
	Dialog.addToSameRow();
	Dialog.addNumber("Resize factor ", marker_scale, 3, 5,"");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Preprocessing:\n", 13, "#f26852");	
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Gaussian Blur Radius", cell_filter_scale, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Segmentation:\n", 13, "#f26852");	
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Segmentation Method", segmentation_methods, cell_method);
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("> Threshold: Method", threshold_methods, cell_threshold);
	Dialog.addToSameRow();
	Dialog.addNumber("> Fixed, threshold value", cell_fixed_threshold_value, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("> Stardist: Probability", stardist_probability, 2, 4, "");
	Dialog.addToSameRow();
	Dialog.addNumber("Overlap", stardist_overlap, 2, 4, "");
	Dialog.addNumber("Prominence", prominence, 0, 4, "");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Min. ROI Area", cell_area_min, 0, 5, micron+"2");
	Dialog.addToSameRow();
	Dialog.addNumber("Max. ROI Area", cell_area_max, 0, 5, micron+"2");
	Dialog.show();

	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
	
	suffix_time						= Dialog.getChoice();		print("Image Type:",suffix_time);
	pixel_size						= Dialog.getNumber(); 		print("Pixel Size:",pixel_size);
	roi_type						= Dialog.getChoice();		print("ROI type:",roi_type);
	marker_identifier				= Dialog.getString();		print("Additional Marker:", marker_identifier);
	suffix_marker					= Dialog.getChoice();		print("Additional Marker:", suffix_marker);
	marker_scale					= Dialog.getNumber();		print("Additional Marker:", marker_scale);
	cell_filter_scale				= Dialog.getNumber(); 		print("Filter Scale:",cell_filter_scale);
	cell_method						= Dialog.getChoice();		print("Segmentation Method:", cell_method);
	cell_threshold					= Dialog.getChoice(); 		print("Threshold:",cell_threshold);
	cell_fixed_threshold_value		= Dialog.getNumber(); 		print("Fixed Threshold Value:",cell_fixed_threshold_value);
	stardist_probability			= Dialog.getNumber();		print("Stardist probability:",stardist_probability);
	stardist_overlap				= Dialog.getNumber();		print("Stardist overlap:",stardist_overlap);
	prominence						= Dialog.getNumber(); 		print("Prominence:",prominence);
	cell_area_min					= Dialog.getNumber(); 		print("Min. ROI Diameter:",cell_area_min);
	cell_area_max					= Dialog.getNumber(); 		print("Max. ROI Diameter:",cell_area_max);
}

function measureStack(id,roi_type, prefix)  
{
	selectImage(id);
	getDimensions(width, height, channels, slices, frames);
	if(roi_type=="Full FOV") // whole image
	{
		run("Select All");
		roiManager("add");
	}
	if(roi_type!="Full FOV")
	{
		if(roi_type=="Additional Marker"){
			path = dir+prefix+marker_identifier+suffix_marker;
			if(File.exists(path)){
				run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
				idm=getImageID();
				run("Scale...", "x="+marker_scale+" y="+marker_scale+" width="+width+" height="+height+" interpolation=Bilinear average create");
				avg_id = getImageID;
				selectImage(idm); close();
				selectImage(avg_id);
				calibrateImage(id);
				rename("avg_id");
			}else {
				print(path + " not found.");
				exit;
			}
		}else{
			run("Z Project...", "start=1 projection=[Average Intensity]");
			avg_id = getImageID;
			selectImage(avg_id);
			rename("avg_id");
		}
		run("Gaussian Blur...", "sigma="+cell_filter_scale);
		if(cell_method=="Stardist"){
			selectImage(avg_id);
			run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], "
			+"args=['input':"+'avg_id'+", 'modelChoice':'Versatile (fluorescent nuclei)',"
			+"'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99',"
			+"'probThresh':'"+stardist_probability+"', 'nmsThresh':'"+stardist_overlap+"', 'outputType':'ROI Manager', 'nTiles':'1', "
			+"'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', "
			+"'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
			
			selectImage(avg_id); close();
			newImage("avg_id", "8-bit black", width, height, 1);
			avg_id = getImageID;
			calibrateImage(avg_id);
			nr_rois = roiManager("count");
			for(r = 0; r < nr_rois; r++)
			{
				roiManager("select",r);
				run("Enlarge...", "enlarge=1");
				run("Clear");
				run("Enlarge...", "enlarge=-1");
				run("Fill");
			}
			roiManager("Deselect");
			roiManager("reset");
			setThreshold(1,255);
			setAutoThreshold("Default dark");
		}
		else if(cell_threshold == "Fixed")
		{
			setAutoThreshold("Default dark");
			getThreshold(min_thresh,max_thresh); 
			setThreshold(cell_fixed_threshold_value,max_thresh);
		}
		else 
		{
			setAutoThreshold(cell_threshold+" dark"); 
		}
	}
	
	if(roi_type=="Calcium Signal" || roi_type=="Additional Marker") // single cells
	{
		selectImage(avg_id);
		run("Select None");
		run("Find Maxima...", "prominence="+prominence+" above output=[Segmented Particles]");
		max_id = getImageID;
		selectImage(max_id);
		run("Analyze Particles...", "size="+cell_area_min+"-"+cell_area_max+" add");
		selectImage(max_id);
		close();
		selectImage(avg_id); 
		close();
	}
	selectImage(id);
	roiManager("Deselect");
	roiManager("Multi Measure");
}

function createOverlay(names)
{
	setForegroundColor(25, 25, 25);
	fields = names.length;
	print(fields,"images");
	for(i=0;i<fields;i++)
	{
		prefix = names[i];
		print(i+1,"/",fields,":",prefix);
		path = dir+prefix+suffix_time;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		idTime=getImageID();
		run("Z Project...", "start=1 projection=[Average Intensity]");
		id = getImageID;
		selectImage(idTime); close();
		selectImage(id);
		Stack.getDimensions(w,h,channels,slices,frames); 
		
		if(!Stack.isHyperStack && channels == 1 && channels*slices*frames!=1)
		{
			channels = slices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
		}

		id = getImageID;
		setSlice(nSlices);
		run("Add Slice","add=channel");
		file_roi=output_dir+prefix+"_roi_set.zip"; 
		print(file_roi);
		if(File.exists(file_roi))
		{	
			selectImage(id);
			setSlice(nSlices);
			roiManager("Open",file_roi);
			roiManager("deselect");
			roiManager("Fill");
			roiManager("reset");
		}
		Stack.getDimensions(w,h,channels,slices,frames); 
		if(!Stack.isHyperStack && channels == 1)
		{
			channels = slices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
		}
	}
	run("Concatenate...", "all_open title=[Concatenated Stacks]");
	Stack.getDimensions(w,h,newchannels,slices,frames);
	for(c=1;c<=channels;c++){Stack.setChannel(c);Stack.setFrame(round(frames/2));resetMinAndMax;}
	range = pow(2,bitDepth);
	for(c=channels+1;c<=newchannels;c++){Stack.setChannel(c);setMinAndMax(0,range/2);}
	Stack.setChannel(newchannels);
	run("Grays");
	run("Make Composite");
	
	if(roi_type=="Additional Marker"){
		for(i=0;i<fields;i++)
		{
			prefix = names[i];
			print(i+1,"/",fields,":",prefix);
			path = dir+prefix+marker_identifier+suffix_marker;
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
			idOrig=getImageID();
			getDimensions(width, height, channels, slices, frames);
			run("Scale...", "x="+marker_scale+" y="+marker_scale+" width="+width+" height="+height+" interpolation=Bilinear average create");
			id = getImageID;
			selectImage(idOrig); close();
			selectImage(id);
			Stack.getDimensions(w,h,channels,slices,frames); 
			
			if(!Stack.isHyperStack && channels == 1 && channels*slices*frames!=1)
			{
				channels = slices;
				run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
			}
	
			id = getImageID;
			setSlice(nSlices);
			run("Add Slice","add=channel");
			file_roi=output_dir+prefix+"_roi_set.zip"; 
			print(file_roi);
			if(File.exists(file_roi))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",file_roi);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
			Stack.getDimensions(w,h,channels,slices,frames); 
			if(!Stack.isHyperStack && channels == 1)
			{
				channels = slices;
				run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
			}
		}
		
		run("Concatenate...", "all_open title=[Concatenated Stacks]");
		Stack.getDimensions(w,h,newchannels,slices,frames);
		for(c=1;c<=channels;c++){Stack.setChannel(c);Stack.setFrame(round(frames/2));resetMinAndMax;}
		range = pow(2,bitDepth);
		for(c=channels+1;c<=newchannels;c++){Stack.setChannel(c);setMinAndMax(0,range/2);}
		Stack.setChannel(newchannels);
		run("Grays");
		run("Make Composite");
	}
	
	if(roi_type=="Additional Marker"){
		idStack=getImageID();
		Stack.getDimensions(w,h,newchannels,newslices,newframes);
		run("Duplicate...", "duplicate frames=1-"+newframes/2);
		selectImage(idStack);
		run("Duplicate...", "duplicate frames="+newframes/2+1+"-"+newframes);
		selectImage(idStack);close();
	}
}

function toggleOverlay()
{	
	run("Select None"); 
	roiManager("deselect");
	roiManager("Show All without labels");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}
function calibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!=micron)run("Properties...", " unit="+micron+" pixel_width="+pixel_size+" pixel_height="+pixel_size);
	else pixel_size = pixelWidth;
}

function decalibrateImage(id)
{
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(unit!="pixel")run("Properties...", " unit=pixel pixel_width=1 pixel_height=1");
}
