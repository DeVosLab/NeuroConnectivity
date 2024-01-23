/*	NeuroConnectivity_Morph_v01.ijm
	*************4**********
	Author: 			Winnok H. De Vos
	Modified by: 		Marlies Verschuuren
	Date Created: 		Februari 20, 2023
	Date Last Modified:	Januari 23, 2024
 	
 	Description:
 	ImageJ/Fiji macro set to quantify neuronal connectivity in images of iPSC derived- and primary neuronal cultures. 
 	This macro is derived from an Acapella® (PerkinElmer) script (https://github.com/VerschuurenM/NeuronalConnectivity). 
 	Cultures are immunochemically labeled for a nuclear marker, dendrite marker and a pre- and postsynaptic marker. 
 	After maximum projection of acquired z-stacks, nuclei are detected using a manually assigned threshold or by applying the trained convolution neural network Stardist. 
 	Neurites are identified using a rough (user-defined threshold) and fine (user-defined threshold after tubeness filtering) segmentation, which is a simplified version of MorphoNeuroNet. 
 	Next, the nuclei mask is subtracted from the neurite mask after which the neurite mask is dilated to obtain a search region in which the pre- and postsynaptic spots are detected. 
 	The spots are first enhanced using a gaussian, Laplacian or multi-scale Laplacian filter with a user-defined kernel size. 
 	Next, a user-defined threshold is applied to segment the spots. 
 	The resolution of the microscope setup does not allow determining the exact location of individual markers within a synapse, but this is not the intention of the assay. 
 	Instead, the lower resolution is exploited to define synapses as those objects that demonstrate an overlap of minimum 1 pixel between the pre- and postsynaptic spots. 
 	In addition to this object-based colocalization, the Pearson correlation of the pre- and postsynaptic channel is calculated in the search region as an intensity-based colocalization metric.

	Citation
	Please cite these papers if you are using this script in your research:
	* Verschuuren M, Verstraelen P, Garcia-Diaz Barriga G, Cilissen I, Coninx E, Verslegers M, Larsen PH, Nuydens R, De Vos WH. High-throughput microscopy exposes a pharmacological window in which dual leucine zipper kinase inhibition preserves neuronal network connectivity. Acta Neuropathol Commun. 2019 Jun 4;7(1):93. doi: 10.1186/s40478-019-0741-3. Erratum in: Acta Neuropathol Commun. 2019 Aug 14;7(1):131. PMID: 31164177; PMCID: PMC6549294.
	* Verstraelen P, Garcia-Diaz Barriga G, Verschuuren M, Asselbergh B, Nuydens R, Larsen PH, Timmermans JP, De Vos WH. Systematic Quantification of Synapses in Primary Neuronal Culture. iScience. 2020 Sep 7;23(9):101542. doi: 10.1016/j.isci.2020.101542. PMID: 33083769; PMCID: PMC7516133.
	Neurite detection is based on:
	* Pani G, De Vos WH, Samari N, de Saint-Georges L, Baatout S, Van Oostveldt P, Benotmane MA. MorphoNeuroNet: an automated method for dense neurite network analysis. Cytometry A. 2014 Feb;85(2):188-99. doi: 10.1002/cyto.a.22408. Epub 2013 Nov 12. PMID: 24222510. 
	If the Stardist algorithm is applied to detect nuclei, please cite: 
	* Uwe Schmidt, Martin Weigert, Coleman Broaddus, and Gene Myers. Cell Detection with Star-convex Polygons. International Conference on Medical Image Computing and Computer-Assisted Intervention (MICCAI), Granada, Spain, September 2018.
	
	
	Change Log
	v1: * Fix Bug Tube Enhancement & Contrast
	
 	_________________________________________________________________
*/

/*
 	***********************

	Variable initiation

	***********************
*/

//	String variables
var dir								= "";										//	directory
var log_path						= "";										//	path for the log file
var micron							= getInfo("micrometer.abbreviation");		// 	micro symbol
var neurites_roi_set 				= "";										//	neurites ROIsets
var neurites_results 				= "";										//	neurites results
var search_roi_set 					= "";										//	neurites search ROIsets
var search_results 					= "";										//	neurites search results
var nuclei_preprocessing_method		= "Median";									// 	method used for preprocessing nuclei image
var nuclei_roi_set 					= "";										//	nuclear ROIsets
var nuclei_results 					= "";										//	nuclear results
var nuclei_segmentation_method		= "Stardist";								// 	Method used for nuclei detection
var nuclei_threshold				= "Fixed";									//	threshold method for segmentation of nuclei 
var order							= "xyczt(default)";							//	hyperstack dimension order
var output_dir						= "";										//	dir for analysis output
var results							= "";										//	summarized results	
var skeleton_results				= "";										//	skeleton results
var skeleton_roi_set				= "";										// 	skeleton roi 
var spots_a_results					= "";										//	spot results ch A
var spots_a_roi_set					= "";										//	spot ROI sets ch A			
var spots_a_segmentation_method		= "Laplace";								//	spots segmentation method ch A
var spots_a_threshold				= "Fixed";									//	threshold method for spot segmentation ch A
var spots_b_results					= "";										//	spot results ch B
var spots_b_roi_set					= "";										//	spot ROI sets ch B	
var spots_b_segmentation_method		= "Laplace";								//	spots segmentation method ch B
var spots_b_threshold				= "Fixed";									//	threshold method for spot segmentation ch B
var suffix							= ".tif";									//	suffix for specifying the file type

//	Number variables
var channels						= 3;										//	number of channels
var fields							= 4;										//	number of xy positions
var image_height					= 1000;										//	image height
var image_width						= 1000;										//	image width
var slices							= 7;										//	number of z-slices
var neurites_channel				= 1;										
var neurites_iteration				= 2;										//  pixel dilation search region
var neurites_enhance_radius			= 2;										// 	radius enhancement method neurite segmentation
var neurites_median_radius			= 2;										//  radius median filter neurites segmentation
var neurites_min_area				= 100;										//  minimal area size neurite fragments
var neurites_nuclei_dilation		= 15; 										// 	pixel dilation nuclei mask	
var neurites_threshold_fine			= 250;										//  threshold fine neurite segmentation
var neurites_threshold_rough		= 4000;										//  threshold rough neurite segmentation
var nuclei_channel					= 3;										//	channel used for segmentation of nuclei 
var nuclei_filter_scale				= 8;										// 	radius for nuclei smoothing/laplace
var nuclei_fixed_threshold_value	= 100;										//	fixed maximum threshold for nuclei segmentation (if auto doesn't work well);
var nuclei_min_circularity			= 0.0;										//	min circularity
var nuclei_min_area					= 15;										//	calibrated min nuclear size (in µm2)
var nuclei_max_area					= 500;										//	calibrated max nuclear size (in µm2)
var nuclei_min_intensity			= 1500;										// 	Min nuc intensity to avoid false positive (stardist)
var nuclei_overlap 					= 0.3;										//	nuclei_overlap amount tolerated for Stardist nuclei detection
var nuclei_probability				= 0.15;										//	minimal nuclei_probability for Stardist nuclei detection
var pixel_size						= 0.18;										//	pixel size (µm)
var spots_a_channel					= 2;										//	channel A used for segmentation of spots
var spots_a_filter_scale			= 2;										//	scale for Laplacian spots ch A
var spots_a_fixed_threshold_value	= 250;										//	fixed maximum threshold for spot segmentation 
var spots_a_max_area				= 50;										//	max spot size in pixels ch A
var spots_a_min_area				= 2;										//	min spot size in pixels ch A
var spots_b_channel					= 0;										//	channel B used for segmentation of spots
var spots_b_filter_scale			= 2;										//	scale for Laplacian spots ch B
var spots_b_fixed_threshold_value	= 400;										//	fixed maximum threshold for spot segmentation
var spots_b_max_area				= 50;										//	max spot size in pixels ch B
var spots_b_min_area				= 2;										//	min spot size in pixels ch B

//	Boolean Parameters
var colocalize_spots				= false;									//	colocalize spot ROIs
var neurites_clahe					= true;										//  local contrast enhancement for neurite segmentation
var neurites_median					= true;										// 	apply median filter prior to neurite segmentation
var nuclei_background				= true;										//	subtract nuclei_background for nuclei segmentation
var nuclei_clahe					= true;										// 	local contrast enhancement for nuclei segmentation
var nuclei_watershed				= false;									//	use nuclei_watershed for nuclei segmentation
var segment_nuclei					= true;										// 	segment nuclei
var segment_neurites				= true;										// 	segment neurites
var segment_spots					= true;										// 	segment spots
var z_project						= true;										//	generation of max projection images

//	Arrays
var dimensions						= newArray("xyczt(default)","xyctz","xytcz","xytzc","xyztc","xyzct");		
var file_list						= newArray(0);								
var file_types 						= newArray(".tif",".tiff",".nd2",".ids",".jpg",".mvd2",".czi");		
var nuclei_segmentation_methods		= newArray("Thresholding","Stardist");
var nuclei_preprocessing_methods	= newArray("Gauss","Median", "Laplace");
var objects 						= newArray("Nuclei","Cells","Spots_a","Spots_b","Folds","Plaques");
var prefixes 						= newArray(0);
var spot_segmentation_methods		= newArray("Gauss","Laplace","Multi-Scale");
var threshold_methods				= getList("threshold.methods");	
var threshold_methods				= Array.concat(threshold_methods,"Fixed");	

/*
 	***********************

		Macros

	***********************
*/

macro "Autorun"
{
	erase(1);
}


macro "Z Project Action Tool - C888 R0099 R3399 R6699 R9999"
{
	setBatchMode(true);
	projectImages();
	setBatchMode("exit and display");
}

macro "Setup Action Tool - C888 T5f16S"
{
	setup();
}

macro "Segment Nuclei Action Tool - C999 V227d"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	idMaskNuc = segmentNuclei(id,nuclei_channel,1); 
	nuclei_nr = roiManager("count");
	if(nuclei_nr<=0){
		selectImage(id);
		getDimensions(w, h, c, s, f);
		newImage("MaskNuc", "8-bit black", w, h, 1);
		idMaskNuc=getImageID();
	}
	selectImage(idMaskNuc); close;
	roiManager("show all without labels");
	setBatchMode("exit and display");
}

macro "Segment Neurites Action Tool - C999 H1157f2c7ffac1f3a1100"
{
	erase(0);
	setBatchMode(true);
	id = getImageID;
	
	if(segment_nuclei){
		idMaskNuc = segmentNuclei(id,nuclei_channel,1); 
		nuclei_nr = roiManager("count");
		if(nuclei_nr<=0){
			selectImage(id);
			getDimensions(w, h, c, s, f);
			newImage("MaskNuc", "8-bit black", w, h, 1);
			idMaskNuc=getImageID();
		}
	}else{
		selectImage(id);
		getDimensions(w, h, c, s, f);
		newImage("MaskNuc", "8-bit black", w, h, 1);
		idMaskNuc=getImageID();
		nuclei_nr=0;
	}
	
	//Detect neurites and create search region. If not create whole FOV selection (with exclusion of nuclei).
	if(segment_neurites){
		roiManager("reset");
		idMaskNeurites = segmentNeurites(id,neurites_channel,idMaskNuc);
		
		//Skeletonize
		roiManager("reset");
		idSkel=skeletonize(idMaskNeurites);
			
		//Search region
		roiManager("reset");
		selectImage(idMaskNeurites);
		run("Duplicate...", " ");
		rename("MaskSearch");
		idMaskSearch=getImageID();
		for (i = 0; i < neurites_iteration; i++) {
			run("Dilate");
		}
		
		//Show ROIs
		selectImage(idMaskNeurites);
		run("Create Selection"); 
		if(selectionType>0){
			roiManager("Add");
			roiManager("Select",RoiManager.size-1);
			roiManager("Set Color", "yellow");
		}	
		selectImage(idMaskNeurites); close;
		
		selectImage(idSkel);
		run("Create Selection"); 
		if(selectionType>0){
			roiManager("Add");
			roiManager("Select",RoiManager.size-1);
			roiManager("Set Color", "cyan");
		}	
		selectImage(idSkel); close;
		
		selectImage(idMaskSearch);
		run("Create Selection"); 
		if(selectionType>0){
			roiManager("Add");
			roiManager("Select",RoiManager.size-1);
			roiManager("Set Color", "magenta");
		}	
	}else{
		roiManager("reset");
		selectImage(id);
		getDimensions(w, h, c, s, f);
		newImage("MaskNeurites", "8-bit white", w, h, 1);
		run("Invert LUT");
		idMaskNeurites=getImageID();
		if(segment_nuclei){
			//Exclude nuclei
			selectWindow("MaskNuc");
			run("Duplicate...", "title=MaskNucDup ");
			for (i = 0; i < neurites_nuclei_dilation; i++) {
				run("Dilate");
			}
			imageCalculator("Subtract create", "MaskNeurites","MaskNucDup");
			rename("MaskSearch");
			idMaskSearch=getImageID();
			run("Create Selection");
			if(selectionType>0){
				roiManager("Add");
				roiManager("Select",0);
				roiManager("Set Color", "#ff00ff");
			}	
			selectWindow("MaskNeurites"); close;
		}else{
			rename("MaskNeurites");
			idMaskSearch=getImageID();
			run("Select All");
			roiManager("Add");
			roiManager("Select",0);
			roiManager("Set Color", "#ff00ff");
		}
	}
	selectImage(idMaskNuc);close;	
	selectImage(idMaskSearch); close;
	roiManager("show all without labels");
	setBatchMode("exit and display");
}

macro "Segment Spots Action Tool - C999 V3633 V4b33 V7633 Va933"
{
	erase(0);
	setBatchMode(true);
	id 		= getImageID;
	cSpot = getNumber("Spot Channel",spots_a_channel);
	
	//Detect nuclei. If not or 0 nuclei detected, create empty mask
	if(segment_nuclei){
		idMaskNuc = segmentNuclei(id,nuclei_channel,1); 
		nuclei_nr = roiManager("count");
		if(nuclei_nr<=0){
			getDimensions(w, h, c, s, f);
			newImage("MaskNuc", "8-bit black", w, h, 1);
			idMaskNuc=getImageID();
		}
	}else{
		getDimensions(w, h, c, s, f);
		newImage("MaskNuc", "8-bit black", w, h, 1);
		idMaskNuc=getImageID();
		nuclei_nr=0;
	}
	//Detect neurites and create search region. If not create whole FOV selection (with exclusion of nuclei).
	
	if(segment_neurites){
		roiManager("reset");
		idMaskNeurites = segmentNeurites(id,neurites_channel,idMaskNuc);
			
		//Search region
		roiManager("reset");
		selectImage(idMaskNeurites);
		run("Duplicate...", " ");
		rename("MaskSearch");
		idMaskSearch=getImageID();
		for (i = 0; i < neurites_iteration; i++) {
			run("Dilate");
		}
		selectImage(idMaskNeurites); close;
	}else{
		roiManager("reset");
		selectImage(id);
		getDimensions(w, h, c, s, f);
		newImage("MaskNeurites", "8-bit white", w, h, 1);
		run("Invert LUT");
		idMaskNeurites=getImageID();
	
		//Search region
		if(nuclei_nr>0){
			selectWindow("MaskNuc");
			run("Duplicate...", "title=MaskNucDup ");
			for (i = 0; i < neurites_nuclei_dilation; i++) {
				run("Dilate");
			}
			imageCalculator("Subtract create", "MaskNeurites","MaskNucDup");
			rename("MaskSearch");
			idMaskSearch=getImageID();
			selectImage(idMaskSearch);
			selectImage(idMaskNeurites); close;
			selectWindow("MaskNucDup"); close;
		}else{
			selectImage(idMaskNeurites);
			idMaskSearch=getImageID();
			rename("MaskSearch");
		}
	}
	
	selectImage(idMaskNuc);close;	
	
	roiManager("reset");
	if(cSpot == spots_a_channel){args = newArray(spots_a_channel,spots_a_segmentation_method,spots_a_threshold,spots_a_fixed_threshold_value,spots_a_filter_scale,spots_a_min_area,spots_a_max_area);}
	if(cSpot== spots_b_channel){args = newArray(spots_b_channel,spots_b_segmentation_method,spots_b_threshold,spots_b_fixed_threshold_value,spots_b_filter_scale,spots_b_min_area,spots_b_max_area);}
	snr = segmentSpots(id,cSpot,idMaskSearch,args);
	
	selectImage(idMaskSearch);
	run("Create Selection"); 
	if(selectionType>0){
		roiManager("Add");
		roiManager("Select",snr);
		roiManager("Set Color", "#ff00ff");
	}
	
	selectImage(idMaskSearch); close;
	selectImage(id);
	roiManager("show all without labels");
	setBatchMode("exit and display");
}

macro "Analyse Single Image Action Tool - C888 T5f161"
{
	erase(0);
	setBatchMode(true);
	dir = getInfo("image.directory");
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir)){
		File.makeDirectory(output_dir);
	}else{
		overwrite=getBoolean("There is already an output folder. Delete previous results?", "Yes", "No");
		if(overwrite){
			list_output = getFileList(output_dir);
			for (i=0; i<list_output.length; i++){
				File.delete(output_dir+list_output[i]);
			}
			File.delete(output_dir);
			File.makeDirectory(output_dir);
		}
	}
	start = getTime();
	title = getTitle; 
	prefix = substring(title,0,lastIndexOf(title,suffix));
	setFileNames(prefix);
	id = getImageID;
	idMaskNeurites=0;
	idMaskSearch=0;
	
	//Detect nuclei. If not or 0 nuclei detected, create empty mask
	if(segment_nuclei){
		idMaskNuc = segmentNuclei(id,nuclei_channel,1); 
		nuclei_nr = roiManager("count");
		if(nuclei_nr>0){
			roiManager("Save",nuclei_roi_set);
		}else{
			getDimensions(w, h, c, s, f);
			newImage("MaskNuc", "8-bit black", w, h, 1);
			idMaskNuc=getImageID();
		}
	}else{
		nuclei_nr=0;
		getDimensions(w, h, c, s, f);
		newImage("MaskNuc", "8-bit black", w, h, 1);
		idMaskNuc=getImageID();
	}
	
	//Detect neurites and create search region. If not create whole FOV selection (with exclusion of nuclei).
	if(segment_neurites){
		roiManager("reset");
		idMaskNeurites = segmentNeurites(id,neurites_channel,idMaskNuc);
		if(RoiManager.size>0){
			roiManager("Save",neurites_roi_set);
			
			//Skeletonize
			roiManager("reset");
			idSkel=skeletonize(idMaskNeurites);
			roiManager("Save",skeleton_roi_set);
			
			//Search region
			roiManager("reset");
			selectImage(idMaskNeurites);
			run("Duplicate...", " ");
			rename("MaskSearch");
			idMaskSearch=getImageID();
			for (i = 0; i < neurites_iteration; i++) {
				run("Dilate");
			}
			run("Create Selection"); 
			if(selectionType>0){
				roiManager("Add");
				roiManager("Save",search_roi_set);
			}
		}else{
			selectImage(id);
			getDimensions(w, h, c, s, f);
			newImage("MaskSearch", "8-bit black", w, h, 1);
			idMaskSearch=getImageID();
		}
		
		selectImage(idMaskNeurites); close;
		selectImage(idSkel); close;
	}else{
		roiManager("reset");
		selectImage(id);
		getDimensions(w, h, c, s, f);
		newImage("MaskNeurites", "8-bit white", w, h, 1);
		run("Invert LUT");
		idMaskNeurites=getImageID();
			
		//Search region
		if(nuclei_nr>0){
			selectWindow("MaskNuc");
			run("Duplicate...", "title=MaskNucDup ");
			for (i = 0; i < neurites_nuclei_dilation; i++) {
				run("Dilate");
			}
			imageCalculator("Subtract create", "MaskNeurites","MaskNucDup");
			rename("MaskSearch");
			idMaskSearch=getImageID();
			selectImage(idMaskSearch);
			run("Create Selection");
			roiManager("Add");
			roiManager("Save",search_roi_set);
			selectImage(idMaskNeurites); close;
			selectWindow("MaskNucDup"); close;
		}else{
			selectImage(idMaskNeurites);
			idMaskSearch=getImageID();
			rename("MaskSearch");
			run("Select All");
			roiManager("Add");
			roiManager("Save",search_roi_set);
		}
	}
	
	selectImage(idMaskNuc);close;
	
	roiManager("reset");
	if(segment_spots && spots_a_channel>0)
	{
		args	= newArray(spots_a_channel,spots_a_segmentation_method,spots_a_threshold,spots_a_fixed_threshold_value,spots_a_filter_scale,spots_a_min_area,spots_a_max_area);
		snr 	= segmentSpots(id,spots_a_channel,idMaskSearch,args);
		if(snr>0)
		{
			roiManager("Save",spots_a_roi_set);
			roiManager("reset");
		}
	}
	if(segment_spots && spots_b_channel>0)
	{
		args	= newArray(spots_b_channel,spots_b_segmentation_method,spots_b_threshold,spots_b_fixed_threshold_value,spots_b_filter_scale,spots_b_min_area,spots_b_max_area);
		snr 	= segmentSpots(id,spots_b_channel,idMaskSearch,args);
		if(snr>0)
		{
			roiManager("Save",spots_b_roi_set);
			roiManager("reset");
		}
	}
	if(isOpen(idMaskSearch)){	
		selectImage(idMaskSearch); close;
	}
	
	readout = analyzeRegions(id);
	
	if(readout)summarizeResults(id);

	run("Select None");
	
	print((getTime()-start)/1000,"sec");
	print("Analysis Done");
	setBatchMode("exit and display");
}

macro "Batch Analysis Action Tool - C888 T5f16#"
{
	erase(1);
	setBatchMode(true);
	setDirectory(true);
	prefixes = scanFiles();
	fields = prefixes.length;
	setup();
	start = getTime();
	for(x=0;x<fields;x++)
	{
		prefix = prefixes[x];
		file = prefix+suffix;
		setFileNames(prefix);
		print(x+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		//open(path);
		id = getImageID;
		idMaskNeurites=0;
		idMaskSearch=0;
	
		//Detect nuclei. If not or 0 nuclei detected, create empty mask
		if(segment_nuclei){
			idMaskNuc = segmentNuclei(id,nuclei_channel,1); 
			nuclei_nr = roiManager("count");
			if(nuclei_nr>0){
				roiManager("Save",nuclei_roi_set);
			}else{
				getDimensions(w, h, c, s, f);
				newImage("MaskNuc", "8-bit black", w, h, 1);
				idMaskNuc=getImageID();
			}
		}else{
			nuclei_nr=0;
			getDimensions(w, h, c, s, f);
			newImage("MaskNuc", "8-bit black", w, h, 1);
			idMaskNuc=getImageID();
		}
		
		//Detect neurites and create search region. If not create whole FOV selection (with exclusion of nuclei).
		if(segment_neurites){
			roiManager("reset");
			idMaskNeurites = segmentNeurites(id,neurites_channel,idMaskNuc);
			if(RoiManager.size>0){
				roiManager("Save",neurites_roi_set);
				
				//Skeletonize
				roiManager("reset");
				idSkel=skeletonize(idMaskNeurites);
				roiManager("Save",skeleton_roi_set);
				selectImage(idSkel); close;
				
				//Search region
				roiManager("reset");
				selectImage(idMaskNeurites);
				run("Duplicate...", " ");
				rename("MaskSearch");
				idMaskSearch=getImageID();
				for (i = 0; i < neurites_iteration; i++) {
					run("Dilate");
				}
				run("Create Selection"); 
				if(selectionType>0){
					roiManager("Add");
					roiManager("Save",search_roi_set);
				}
			}else{
				selectImage(id);
				getDimensions(w, h, c, s, f);
				newImage("MaskSearch", "8-bit black", w, h, 1);
				idMaskSearch=getImageID();
			}
			selectImage(idMaskNeurites); close;

		}else{
			roiManager("reset");
			selectImage(id);
			getDimensions(w, h, c, s, f);
			newImage("MaskNeurites", "8-bit white", w, h, 1);
			run("Invert LUT");
			idMaskNeurites=getImageID();
				
			//Search region
			if(nuclei_nr>0){
				selectWindow("MaskNuc");
				run("Duplicate...", "title=MaskNucDup ");
				for (i = 0; i < neurites_nuclei_dilation; i++) {
					run("Dilate");
				}
				imageCalculator("Subtract create", "MaskNeurites","MaskNucDup");
				rename("MaskSearch");
				idMaskSearch=getImageID();
				selectImage(idMaskSearch);
				run("Create Selection");
				roiManager("Add");
				roiManager("Save",search_roi_set);
				selectImage(idMaskNeurites); close;
				selectWindow("MaskNucDup"); close;
			}else{
				selectImage(idMaskNeurites);
				idMaskSearch=getImageID();
				rename("MaskSearch");
				run("Select All");
				roiManager("Add");
				roiManager("Save",search_roi_set);
			}
		}
		
		selectImage(idMaskNuc);close;
		
		roiManager("reset");
		if(segment_spots && spots_a_channel>0)
		{
			args	= newArray(spots_a_channel,spots_a_segmentation_method,spots_a_threshold,spots_a_fixed_threshold_value,spots_a_filter_scale,spots_a_min_area,spots_a_max_area);
			snr 	= segmentSpots(id,spots_a_channel,idMaskSearch,args);
			if(snr>0)
			{
				roiManager("Save",spots_a_roi_set);
				roiManager("reset");
			}
		}
		if(segment_spots && spots_b_channel>0)
		{
			args	= newArray(spots_b_channel,spots_b_segmentation_method,spots_b_threshold,spots_b_fixed_threshold_value,spots_b_filter_scale,spots_b_min_area,spots_b_max_area);
			snr 	= segmentSpots(id,spots_b_channel,idMaskSearch,args);
			if(snr>0)
			{
				roiManager("Save",spots_b_roi_set);
				roiManager("reset");
			}
		}
		if(isOpen(idMaskSearch)){	
			selectImage(idMaskSearch); close;
		}
		
		readout = analyzeRegions(id);
		
		if(readout)summarizeResults(id);
		
		run("Select None");
		selectImage(id); close;
		erase(0);
	}
	print((getTime()-start)/1000,"sec");
	if(isOpen("Log")){selectWindow("Log");saveAs("txt",log_path);}
	concatenateResults();
	saveAs(".txt",output_dir+"ConcatenatedResults.txt");
	print("Complete Analysis Done");
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
	setDirectory(false);
	prefixes = scanFiles();
	subset = getBoolean(""+prefixes.length+" images, select subset?");
	if(subset)
	{
		for (i = 0; i < prefixes.length; i++) {
			print(i+1,"/",prefixes.length,":",prefixes[i]);
		}
		Dialog.create("Verification Stack");
		Dialog.addMessage("Select wells to visualize...")
		Dialog.addNumber("Start file index", 1);
		Dialog.addNumber("Stop file index", prefixes.length);
  		Dialog.show;
  		
		startIndex=Dialog.getNumber()-1;
		stopIndex=Dialog.getNumber();
		names =  Array.slice(prefixes, startIndex, stopIndex); 
	}
	else names = prefixes;
	createOverlay(names);
	setBatchMode("exit and display");
	run("Tile");
	run("Channels Tool... ");
}

/*
 	***********************

		Functions

	***********************
*/

function projectImages()
{
	erase(1);
	Dialog.create("Project Images...");
	Dialog.addString("Destination Directory Name","Images",25);
	Dialog.addString("Add a prefix","Image",25);
	Dialog.addChoice("Import format",file_types,".nd2");
	Dialog.addChoice("Export format",file_types,suffix);
	Dialog.show;
	dest 		= Dialog.getString;
	pre			= Dialog.getString;
	ext			= Dialog.getChoice;
	suffix 		= Dialog.getChoice;
	dir 		= getDirectory("");
	file_list 	= getFileList(dir);
	destination_dir 	= dir+dest+File.separator;
	File.makeDirectory(destination_dir);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,ext))
		{		
			print(i+1);
			run("Bio-Formats Importer", "open=["+path+"] color_mode=Default concatenate_series open_all_series view=Hyperstack ");
			ns = nSlices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+ns/channels+" frames=1 display=Color");
			id = getImageID;
			title = getTitle;
			run("Z Project...", "projection=[Max Intensity]");
			zid = getImageID;		
			selectImage(zid); saveAs(suffix,destination_dir+pre+title+suffix);
			selectImage(id); close;
			selectImage(zid); close;
		}
	}
	print("Done");
}

function setOptions()
{
	run("Options...", "iterations=1 count=1");
	run("Colors...", "foreground=white nuclei_background=black selection=yellow");
	run("Overlay Options...", "stroke=red width=1 fill=none");
	setBackgroundColor(0, 0, 0);
	setForegroundColor(255,255,255);
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
	run("Select None"); 
	run("Remove Overlay");
	if(all){
		print("\\Clear");
		run("Close All");
	}
	run("Clear Results");
	roiManager("reset");
	run("Collect Garbage");
}

function setDirectory(delete)
{
	dir = getDirectory("Choose a Source Directory");
	file_list = getFileList(dir);
	output_dir = dir+"Output"+File.separator;
	if(!File.exists(output_dir)){
		File.makeDirectory(output_dir);
	}else if (delete){
		overwrite=getBoolean("There is already an output folder. Delete previous results?", "Yes", "No");
		if(overwrite){
			list_output = getFileList(output_dir);
			for (i=0; i<list_output.length; i++){
				File.delete(output_dir+list_output[i]);
			}
			File.delete(output_dir);
			File.makeDirectory(output_dir);
		}
	}
	log_path = output_dir+"Log.txt";
}

function setFileNames(prefix)
{
	nuclei_roi_set 				= output_dir+prefix+"_nuclei_roi_set.zip";
	nuclei_results 				= output_dir+prefix+"_nuclei_results.txt";
	neurites_roi_set 			= output_dir+prefix+"_neurites_roi_set.zip";
	neurites_results			= output_dir+prefix+"_neurites_results.txt";
	search_roi_set 				= output_dir+prefix+"_neurites_search_roi_set.zip";
	search_results				= output_dir+prefix+"_neurites_search_results.txt";
	skeleton_roi_set			= output_dir+prefix+"_neurites_skeleton_roi_set.zip";
	skeleton_results			= output_dir+prefix+"_neurites_skeleton_results.txt";
	spots_a_roi_set				= output_dir+prefix+"_spots_a_roi_set.zip";
	spots_b_roi_set				= output_dir+prefix+"_spots_b_roi_set.zip";
	spots_a_results				= output_dir+prefix+"_spots_a_results.txt";
	spots_b_results				= output_dir+prefix+"_spots_b_results.txt";
	results						= output_dir+prefix+"_summary.txt";
}

function scanFiles()
{
	prefixes = newArray(0);
	for(i=0;i<file_list.length;i++)
	{
		path = dir+file_list[i];
		if(endsWith(path,suffix) && indexOf(path,"flatfield")<0)
		{
			print(path);
			prefixes = Array.concat(prefixes,substring(file_list[i],0,lastIndexOf(file_list[i],suffix)));			
		}
	}
	return prefixes;
}

function setup()
{
	setOptions();
	Dialog.createNonBlocking("NeuroConnectivity Morph Settings");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("-------------------   General parameters  -------------------", 13, "#f26852");	
	Dialog.setInsets(0,0,0);
	Dialog.addChoice("Image Type", file_types, suffix);
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Pixel Size", pixel_size, 3, 5, micron+"");
	Dialog.setInsets(0,0,0);
	Dialog.addNumber("Channel Number", channels, 0, 1, "");
	Dialog.setInsets(0,0,0);
	Dialog.addMessage("Detectable Objects:");
	labels = newArray(3);defaults = newArray(3);
	labels[0] = "Nuclei";		defaults[0] = segment_nuclei;
	labels[1] = "Neurites";		defaults[1] = segment_neurites;
	labels[2] = "Spots";		defaults[2] = segment_spots;
	Dialog.setInsets(0,110,0);
	Dialog.addCheckboxGroup(1,5,labels,defaults);
	Dialog.show();
	
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
	suffix							= Dialog.getChoice();		print("Image Type:",suffix);
	pixel_size						= Dialog.getNumber(); 		print("Pixel Size:",pixel_size);
	channels 						= Dialog.getNumber();		print("Channels:",channels);
	segment_nuclei					= Dialog.getCheckbox(); 	print("Segment Nuclei:",segment_nuclei);
	segment_neurites				= Dialog.getCheckbox();		print("Segment Cells:",segment_neurites);
	segment_spots 					= Dialog.getCheckbox();		print("Segment Spots:",segment_spots);

	if(segment_nuclei){
		Dialog.createNonBlocking("NeuroConnectivity Morph Settings");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("---------------------------  Nuclei parameters  ---------------------------", 13, "#f26852");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Nuclei are segmented by thresholding, or using a trained classifier (Stardist).\nPreprossing steps and object filtering can be tuned.\n", 12, "#999999");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("General:\n", 13, "#f26852");	
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Nuclei Channel", nuclei_channel, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Segmentation Method", nuclei_segmentation_methods, nuclei_segmentation_method);
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Preprocessing:\n", 13, "#f26852");	
		Dialog.setInsets(0,140,0);	
		Dialog.addCheckbox("Background Subtraction", nuclei_background);
		Dialog.setInsets(0,140,0);	
		Dialog.addCheckbox("Local Contrast Enhancement", nuclei_clahe);
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Preprocessing Method", nuclei_preprocessing_methods, nuclei_preprocessing_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Filter Size", nuclei_filter_scale, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Threshold Settings:\n", 13, "#f26852");	
		Dialog.setInsets(0,0,0);	
		Dialog.addChoice("Threshold Method", threshold_methods, nuclei_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", nuclei_fixed_threshold_value, 0, 4, "");
		Dialog.setInsets(0,140,0);
		Dialog.addCheckbox("Watershed Separation", nuclei_watershed);
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Stardist settings:\n", 13, "#f26852");	
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Probability", nuclei_probability, 2, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Tolerated overlap", nuclei_overlap, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. intensity", nuclei_min_intensity, 2, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Object filters:\n", 13, "#f26852");	
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. Circularity", nuclei_min_circularity, 2, 5, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. Area", nuclei_min_area, 0, 5, micron+"2");
		Dialog.addToSameRow();
		Dialog.addNumber("Max. Area", nuclei_max_area, 0, 5, micron+"2");	
		Dialog.setInsets(20,0,0);
		Dialog.show();
		
		nuclei_channel 					= Dialog.getNumber();		print("Nuclear Channel:",nuclei_channel);
		nuclei_segmentation_method		= Dialog.getChoice();		print("Nuclei Segmentation Method:",nuclei_segmentation_method);
		nuclei_background				= Dialog.getCheckbox();		print("Nuclei Background Subtraction:",nuclei_background);
		nuclei_clahe					= Dialog.getCheckbox();		print("Nuclei Clahe:",nuclei_clahe);
		nuclei_preprocessing_method		= Dialog.getChoice();		print("Nuclei Preprocessing Method:",nuclei_preprocessing_method);
		nuclei_filter_scale				= Dialog.getNumber();		print("Nuclei Filter Scale:",nuclei_filter_scale);
		nuclei_threshold				= Dialog.getChoice();		print("Nuclear Autothreshold:",nuclei_threshold);
		nuclei_fixed_threshold_value	= Dialog.getNumber();		print("Nuclei Fixed Threshold Value:",nuclei_fixed_threshold_value);
		nuclei_watershed 				= Dialog.getCheckbox();		print("Nuclei Watershed:",nuclei_watershed);
		nuclei_probability 				= Dialog.getNumber();		print("Nuclei Stardist - Probability:",nuclei_probability);
		nuclei_overlap 					= Dialog.getNumber();		print("Nuclei Stardist - Overlap:",nuclei_overlap);
		nuclei_min_intensity			= Dialog.getNumber();		print("Nuclei Stardist - Min Intensity:", nuclei_min_intensity);
		nuclei_min_circularity			= Dialog.getNumber();		print("Nuclei Min Circ:",nuclei_min_circularity);
		nuclei_min_area					= Dialog.getNumber();		print("Nuclei Min Nuclear Size:",nuclei_min_area);
		nuclei_max_area					= Dialog.getNumber();		print("Nuclei Max Nuclear Size:",nuclei_max_area);
	}

	if(segment_neurites){
		Dialog.createNonBlocking("NeuroConnectivity Morph Settings");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------  Neurites parameters  -----------------------", 13, "#f26852");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Neurites are detected by a combined rough and fine segmentation. \n", 12, "#999999");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("General :\n", 13, "#f26852");   
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Neurites Channel", neurites_channel,0,4," ");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-----------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Preprocessing:\n", 13, "#f26852");	
		Dialog.setInsets(0,140,0);
		Dialog.addCheckbox("Contrast Enhancement", neurites_clahe);
		Dialog.setInsets(0,140,0);
		Dialog.addCheckbox("Median Filter", neurites_median);
		Dialog.addToSameRow();
		Dialog.addNumber("Median Filter scale ",neurites_median_radius,0,3,"");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Tubeness Enhancement Settings :\n", 13, "#f26852");
		Dialog.setInsets(0,0,0);	
		Dialog.addNumber("Scale ",neurites_enhance_radius,1,3,"");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Threshold Settings :\n", 13, "#f26852");  
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Fixed Threshold Fine", neurites_threshold_fine, 0, 4, "");
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold Rough", neurites_threshold_rough, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Object filters:\n", 13, "#f26852");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min Area Fragments", neurites_min_area, 0, 4, micron+"2");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Search area:\n", 13, "#f26852");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Dilation search region", neurites_iteration, 0, 4, "px");	
		Dialog.addToSameRow();
		Dialog.addNumber("Dilation nuclei exclusion", neurites_nuclei_dilation, 0, 4, "px");	
		Dialog.show();


		
		neurites_channel 				= Dialog.getNumber();		print("Neurites Channel:",neurites_channel);
		neurites_clahe					= Dialog.getCheckbox();		print("Neurites Clahe:",neurites_clahe);
		neurites_median					= Dialog.getCheckbox();		print("Neurites Median:", neurites_median);
		neurites_median_radius			= Dialog.getNumber();		print("Neurites Median radius:", neurites_median_radius);
		neurites_enhance_radius			= Dialog.getNumber();		print("Neurites enhance radius:", neurites_enhance_radius);
		neurites_threshold_fine			= Dialog.getNumber();		print("Neurites threshold fine:", neurites_threshold_fine);
		neurites_threshold_rough		= Dialog.getNumber();		print("Neurites threshold rough:", neurites_threshold_rough);
		neurites_min_area				= Dialog.getNumber();		print("Neurites min area:", neurites_min_area);
		neurites_iteration				= Dialog.getNumber();		print("Neurites dilation search region:", neurites_iteration);
		neurites_nuclei_dilation		= Dialog.getNumber();		print("Neurites dilation nuclei:", neurites_nuclei_dilation);
	}
	
	if(segment_spots){
		Dialog.createNonBlocking("NeuroConnectivity Morph Settings");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("-------------------------  Spot parameters  -------------------------", 13, "#f26852");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Detection of spots in 1 or 2 channels by thresholding after Laplace/Gauss enhancement\nIf there is only one spot channel, set the second to 0.\nIf both are present, colocalize allows detecting reciprocal overlap.", 12, "#999999");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);	
		Dialog.addMessage("Settings Spots A :\n", 13, "#f26852"); 
		Dialog.setInsets(0,0,0);	
		Dialog.addNumber("Spot Channel A", spots_a_channel, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Segmentation Method", spot_segmentation_methods, spots_a_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Filter Scale", spots_a_filter_scale, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Method", threshold_methods, spots_a_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_a_fixed_threshold_value, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. Spot Size", spots_a_min_area, 0, 4, "pixels");
		Dialog.addToSameRow();
		Dialog.addNumber("Max. Spot Size", spots_a_max_area, 0, 4, "pixels");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Settings Spots B :\n", 13, "#f26852"); 
		Dialog.setInsets(0,0,0);	
		Dialog.addNumber("Spot Channel B", spots_b_channel, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Segmentation Method", spot_segmentation_methods, spots_b_segmentation_method);
		Dialog.addToSameRow();
		Dialog.addNumber("Filter Scale", spots_b_filter_scale, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addChoice("Threshold Method", threshold_methods, spots_b_threshold);
		Dialog.addToSameRow();
		Dialog.addNumber("Fixed Threshold", spots_b_fixed_threshold_value, 0, 4, "");
		Dialog.setInsets(0,0,0);
		Dialog.addNumber("Min. Spot Size", spots_b_min_area, 0, 4, "px");
		Dialog.addToSameRow();
		Dialog.addNumber("Max. Spot Size", spots_b_max_area, 0, 4, "px");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("------------------------------------------------------------------", 13, "#dddddd");
		Dialog.setInsets(0,0,0);
		Dialog.addMessage("Colocalisation :\n", 13, "#f26852");	 
		Dialog.setInsets(0,140,0);	
		Dialog.addCheckbox("Colocalize Spots",colocalize_spots);
		Dialog.show();
	
		spots_a_channel					= Dialog.getNumber();		print("Spot A Channel:",spots_a_channel);
		spots_a_segmentation_method		= Dialog.getChoice();		print("Spot A Segmentation Method:",spots_a_segmentation_method);
		spots_a_filter_scale 			= Dialog.getNumber();		print("Spot A Filter Size:",spots_a_filter_scale);
		spots_a_threshold 				= Dialog.getChoice();		print("Spot A AutoThreshold:",spots_a_threshold);
		spots_a_fixed_threshold_value  	= Dialog.getNumber();		print("Spot A Fixed Threshold Value:",spots_a_fixed_threshold_value);
		spots_a_min_area	 			= Dialog.getNumber();		print("Spot A Min. size:", spots_a_min_area);
		spots_a_max_area 				= Dialog.getNumber();		print("Spot A Max. size:",spots_a_max_area);
		
		spots_b_channel					= Dialog.getNumber();		print("Spot B Channel:",spots_b_channel);
		spots_b_segmentation_method		= Dialog.getChoice();		print("Spot B Segmentation Method:",spots_b_segmentation_method);
		spots_b_filter_scale 			= Dialog.getNumber();		print("Spot B Filter Scale:",spots_b_filter_scale);
		spots_b_threshold 				= Dialog.getChoice();		print("Spot B AutoThreshold:",spots_b_threshold);
		spots_b_fixed_threshold_value  	= Dialog.getNumber();		print("Spot B Fixed Threshold Value:",spots_b_fixed_threshold_value);
		spots_b_min_area	 			= Dialog.getNumber();		print("Spot B Min. Size:", spots_b_min_area);
		spots_b_max_area 				= Dialog.getNumber();		print("Spot B Max. Size:",spots_b_max_area);
		
		colocalize_spots 				= Dialog.getCheckbox();		print("Colocalize spot channels",colocalize_spots); 
		
		if(spots_a_channel ==0 || spots_b_channel==0){
			colocalize_spots=false;
		}
	}
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
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

function segmentNuclei(id,c,sel)
{
	// input = multichannel image, output = roiset of nuclear ROIs and if(sel==1) mask incl. border objects
	// output = an image (mid) that contains all ROIs (also touching borders) and roiset of full nuclei
	mid = 0;
	selectImage(id);
	image_width = getWidth;
	image_height = getHeight;
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	cid = getImageID; // the nuclear channel image to be turned into a binary image
	calibrateImage(cid);
	
	selectImage(id);
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	origid = getImageID; 
	
	// preprocess the image
	selectImage(cid);
	if(nuclei_clahe)run("Enhance Local Contrast (CLAHE)", "blocksize=100 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	if(nuclei_background)run("Subtract Background...", "rolling="+round(30/pixel_size));
	if(nuclei_preprocessing_method == "Gauss" && nuclei_filter_scale > 0)run("Gaussian Blur...", "sigma="+nuclei_filter_scale);
	else if(nuclei_preprocessing_method == "Median" && nuclei_filter_scale > 0)run("Median...", "radius="+nuclei_filter_scale);
	else if(nuclei_preprocessing_method == "Laplace")
	{
		run("FeatureJ Laplacian", "compute smoothing="+nuclei_filter_scale); // scale to be adapted depending on nuclei size and SNR
		selectImage(cid); close;
		selectImage("copy Laplacian");
		rename("copy");
		cid = getImageID;
		selectImage(cid);
	}
	if(nuclei_segmentation_method != "Stardist")
	{
		if(nuclei_threshold=="Fixed")
		{
			if(nuclei_preprocessing_method == "Laplace")setAutoThreshold("Default ");
			else setAutoThreshold("Default dark");
			getThreshold(mit,mat); 
			setThreshold(nuclei_fixed_threshold_value,mat);
		}
		else {
			if(nuclei_preprocessing_method == "Laplace")setAutoThreshold(nuclei_threshold); 
			else setAutoThreshold(nuclei_threshold+" dark"); 
		}
		getThreshold(mit,mat); print("Nuclei Threshold:",mit,mat);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		run("Fill Holes");
		if(nuclei_watershed)run("Watershed");
	}
	else if(nuclei_segmentation_method == "Stardist")
	{
		selectImage(cid);
		run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], "
		+"args=['input':"+'copy'+", 'modelChoice':'Versatile (fluorescent nuclei)',"
		+"'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99',"
		+"'probThresh':'"+nuclei_probability+"', 'nmsThresh':'"+nuclei_overlap+"', 'outputType':'ROI Manager', 'nTiles':'1', "
		+"'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', "
		+"'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
		selectImage(cid); close;
		newImage("copy", "8-bit black", image_width, image_height, 1);
		calibrateImage(cid);
		cid = getImageID;
		selectImage(cid);
		selectImage("copy");
		nr_rois = roiManager("count");
		run("Set Measurements...", "mean redirect=None decimal=2");
		for(r = 0; r < nr_rois; r++){
			selectImage(origid);
			roiManager("select",r);
			roiManager("measure");
			if(getResult("Mean", 0)>nuclei_min_intensity){
				selectImage(cid);
				roiManager("select",r);
				run("Enlarge...", "enlarge=1");
				run("Clear");
				run("Enlarge...", "enlarge=-1");
				run("Fill");
			}else{
				roiManager("select",r);
				roiManager("delete");
				r=r-1; 
				nr_rois=nr_rois-1;
			}
			run("Clear Results");
		}
		selectImage(cid);
		roiManager("Deselect");
		roiManager("reset");
		setThreshold(1,255);
		run("Convert to Mask");
	}
	
	run("Set Measurements...", "area centroid mean integrated redirect=None decimal=2");
	if(sel)
	{
		selectImage(cid);
		run("Analyze Particles...", "size="+nuclei_min_area+"-"+nuclei_max_area+" circularity="+nuclei_min_circularity+"-1.00 show=Masks clear include");
		if(isOpen("Mask of copy"))
		{
			selectWindow("Mask of copy"); 
			rename("MaskNuc");
			mid = getImageID;   // the full mask of all particles (for more accurate cell segmentation)
		}
	}
	selectImage(cid);

	run("Analyze Particles...", "size="+nuclei_min_area+"-"+nuclei_max_area+" circularity="+nuclei_min_circularity+"-1.00 show=Nothing exclude clear include add");
	rmc = roiManager("count"); 
	print(rmc,"Nuc. ROI");
	if(rmc==0 && isOpen(mid)){selectImage(mid); close; mid=0;}
	for(i=0;i<rmc;i++)
	{
		roiManager("select",i);
		if(i<9)roiManager("Rename","000"+i+1);
		else if(i<99)roiManager("Rename","00"+i+1);	
		else if(i<999)roiManager("Rename","0"+i+1);	
		else roiManager("Rename",i+1);
	}
	run("Clear Results");
	roiManager("deselect"); 
	roiManager("Measure");
	selectImage(cid); close; 
	selectImage(origid); close;
	return mid;
}

function segmentNeurites(id,c,mid)
{
	selectImage(id);
	run("Select None");
	setSlice(c); 
	run("Duplicate...", "title=copy ");
	cid = getImageID;
	selectImage(cid);
	calibrateImage(cid);
	if(bitDepth==32){
		resetMinAndMax;
		run("16-bit");
	}else{
		resetMinAndMax;
	}
	if(neurites_clahe>0)run("Enhance Local Contrast (CLAHE)", "blocksize=100 histogram=256 maximum=3 mask=*None* fast_(less_accurate)");
	if(neurites_median>0) run("Median...","radius="+neurites_median_radius);
	
	//	Segment fine features
	selectImage(cid);
	run("Tubeness", "sigma="+neurites_enhance_radius+"");
	eid = getImageID; 
	selectImage(eid);
	setAutoThreshold("Default dark");
	getThreshold(mit,mat); 
	print("Original threshold settings neurite fine:",mit,mat);
	setAutoThreshold("Default dark");
	setThreshold(neurites_threshold_fine,maxOf(mat,neurites_threshold_fine));
	run("Convert to Mask");
	rename("Fine");
	

	//	Segment rough features
	selectImage(cid);
	setAutoThreshold("Default dark ");
	getThreshold(mit,mat); 
	print("Original threshold settings neurite rough:",mit,mat);
	setThreshold(neurites_threshold_rough,maxOf(mat,neurites_threshold_rough));	
	run("Convert to Mask");
	rename("Rough");

	//	Combine Fine and Rough
	imageCalculator("OR create", "Rough", "Fine");
	rename("MaskNeurites");
	
	if(segment_nuclei){
		//Exclude nuclei
		selectWindow("MaskNuc");
		run("Duplicate...", "title=MaskNucDup ");
		for (i = 0; i <  neurites_nuclei_dilation; i++) {
			run("Dilate");
		}
		imageCalculator("Subtract create", "MaskNeurites","MaskNucDup");
		rename("MaskNeuritesExclNuc");
		combid=getImageID();
		
		selectWindow("MaskNeurites"); close;
		selectWindow("MaskNucDup"); close;
		selectImage(eid); close;
		selectImage(cid); close;
	}else{
		selectWindow("MaskNeurites");
		combid=getImageID();
		selectImage(eid); close;
		selectImage(cid); close;
	}

	//	Analyze combined mask
	selectImage(combid);
	rename("Combined");
	run("Analyze Particles...", "size="+neurites_min_area+"-Infinity pixel circularity=0-1.00 show=Masks");	
	selectImage(combid); close;
	selectWindow("Mask of Combined");
	rename("MaskNeurites");
	idMaskNeurites = getImageID;
	selectImage(idMaskNeurites);
	run("Create Selection"); 
	if(selectionType>0){
		roiManager("Add");
	}
	
	run("Remove Overlay");
	run("Select None");
	return idMaskNeurites;
}

function skeletonize(idMaskNeurites)
{
	selectImage(idMaskNeurites);
	im_width = getWidth; 
	im_height = getHeight;
	newImage("Skeleton", "8-bit black", im_width, im_height, 1);
	idSkel = getImageID;
	selectImage(idMaskNeurites);
	run("Create Selection"); 
	selectImage(idSkel);
	run("Restore Selection");
	run("Fill");
	run("Select None");
	selectImage(idMaskNeurites);
	run("Select None");
	selectImage(idSkel);
	setThreshold(1,255);
	run("Convert to Mask");
	run("Skeletonize");
	selectImage(idSkel);
	setThreshold(1,255);
	run("Convert to Mask");
	run("Create Selection");
	roiManager("Add");
	return idSkel;
}

function segmentSpots(id,c, idMaskSearch, args){
	spot_channel				= args[0];
	spot_segmentation_method	= args[1];
	spot_threshold_method		= args[2];
	spot_fixed_threshold_value	= args[3];
	scale						= args[4];
	spot_min_area				= args[5];
	spot_max_area				= args[6];
	
	selectImage(id);
	run("Select None");
	if(Stack.isHyperstack)run("Duplicate...", "title=copy duplicate channels="+c);	
	else{setSlice(c);run("Duplicate...","title=copy ");}
	cid = getImageID;
	decalibrateImage(cid);
	selectImage(cid);
	title = getTitle;
	if(spot_segmentation_method=="Gauss")
	{
		run("Duplicate...","title=Gauss");
		run("Gaussian Blur...", "sigma="+scale);
		lap = getImageID;
	}
	else if(spot_segmentation_method=="Laplace")
	{
		run("FeatureJ Laplacian", "compute smoothing="+scale);
		lap = getImageID;
	}
	else if(spot_segmentation_method=="Multi-Scale") 
	{
		e = 0;
		while(e<scale)
		{			
			e++;
			selectImage(cid);
			run("FeatureJ Laplacian", "compute smoothing="+e);
			selectWindow(title+" Laplacian");
			run("Multiply...","value="+e*e);
			rename("scale "+e);
			eid = getImageID;
			if(e>1)
			{
				selectImage(eid);run("Select All");run("Copy");close;
				selectImage(fid);run("Add Slice");run("Paste");
			}
			else fid = getImageID;
		}
		selectImage(fid);
		run("Z Project...", "start=1 projection=[Sum Slices]");
		lap = getImageID;
		selectImage(fid); close;
	}
	selectImage(lap);	
	
	if(spot_threshold_method=="Fixed")
	{
		if(spot_segmentation_method=="Gauss")setAutoThreshold("Default dark");
		else setAutoThreshold("Default ");
		getThreshold(mit,mat); 
		print("Original threshold settings spot SC"+c+":",mit,mat);
		if(spot_segmentation_method=="Gauss")setThreshold(spot_fixed_threshold_value,maxOf(mat,spot_fixed_threshold_value));
		else setThreshold(minOf(mit,-spot_fixed_threshold_value),-spot_fixed_threshold_value);
	}
	else 
	{
		if(spot_segmentation_method=="Gauss")setAutoThreshold(spot_threshold_method + " dark");
		else setAutoThreshold(spot_threshold_method+" ");
	}
	
	run("Analyze Particles...", "size="+spot_min_area+"-"+spot_max_area+" circularity=0.00-1.00 show=Nothing clear include add");
	run("Clear Results");

	for(i=0;i<RoiManager.size;i++)
	{
		selectImage(idMaskSearch);
		roiManager("select", i);
		getStatistics(area, mean, min, max, std, histogram);
		if(mean==0){
			roiManager("delete");
			i=i-1;
		}
	}
	snr = roiManager("count"); print(snr,"spots SC"+c);
	if(snr>10000){print("excessive number, hence reset"); snr=0; roiManager("reset");}
	for(i=0;i<RoiManager.size;i++)
	{
		roiManager("select",i);
		if(i<9)roiManager("Rename","000"+i+1);
		else if(i<99)roiManager("Rename","00"+i+1);	
		else if(i<999)roiManager("Rename","0"+i+1);	
		else roiManager("Rename",i+1);
	}
	//to avoid excessive spot finding when there are no true spots
	selectImage(lap); close;
	selectImage(cid); close;
	selectImage(idMaskSearch); 
	run("Select None");
	run("Remove Overlay");
	return snr;
}

function analyzeRegions(id)
{
	erase(0); 
	mask = 0;
	readout = 0;
	//	analyze neurites rois
	selectImage(id);
	calibrateImage(id);
	
	if(File.exists(search_roi_set)){
		run("Set Measurements...", "area mean standard median min redirect=None decimal=4");
		roiManager("Open",search_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}

		sortResults(); // organize results per channel
		
		saveAs("Measurements",search_results);
		erase(0);
		readout = 1;
	}
	if(segment_neurites && File.exists(neurites_roi_set))
	{
		run("Set Measurements...", "area mean standard median min redirect=None decimal=4");
		roiManager("Open",neurites_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}

		sortResults(); // organize results per channel
		
		saveAs("Measurements",neurites_results);
		erase(0);
		readout = 1;
	}	
	
	if(segment_neurites && File.exists(skeleton_roi_set))
	{
		//ADAPT SKELETON ANALYSE
		selectImage(id);
		w = getWidth;
		h = getHeight;
		newImage("Skeleton", "8-bit Black",w, h, 1);
		idSkel=getImageID();
		roiManager("Open",skeleton_roi_set);
		roiManager("fill");
		setForegroundColor(255, 255, 255);
		run("Convert to Mask");
		run("Analyze Skeleton (2D/3D)", "prune=none");
		selectWindow("Tagged skeleton");
		rename("Tagged");
		idTag = getImageID;
	
		// analyze skeleton elements
		roiManager("reset");
		run("Set Measurements...", "integrated redirect=None decimal=4");
		branch_tot_length 	= 0; 
		branch_av_length	= 0; 
		branch_nr			= 0; 
		junct_nr			= 0; 
		end_nr				= 0; 

		run("Clear Results");
		selectImage(idSkel);
		run("Create Selection"); 
		roiManager("Add");

		// neurites = 127
		selectImage(idTag);
		run("Duplicate...","title=temp");
		idTagDup = getImageID;
		selectImage(idTagDup);
		setThreshold(127,127);
		run("Convert to Mask");
		roiManager("Measure");
		branch_tot_length = getResult("IntDen",0)/255*pixel_size;
		run("Clear Results");
		
		// number branches
		selectImage(idTagDup);
		run("Ultimate Points"); 
		roiManager("Measure");
		branch_nr = getResult("IntDen",0);
		branch_av_length = branch_tot_length/branch_nr;
		run("Clear Results");
		selectImage(idTagDup); close;
		
		// branch points = 70
		selectImage(idTag);
		run("Duplicate...","title=temp");
		idTagDup = getImageID;
		selectImage(idTagDup);
		setThreshold(70,70);
		run("Convert to Mask");
		run("Ultimate Points"); 
		roiManager("Measure");
		junct_nr = getResult("IntDen",0);
		run("Clear Results");
		selectImage(idTagDup); close;
		
		// end points = 30
		selectImage(idTag);
		run("Duplicate...","title=temp");
		idTagDup = getImageID;
		selectImage(idTagDup);
		setThreshold(30,30);
		run("Convert to Mask");
		run("Ultimate Points"); 
		roiManager("Measure");
		end_nr = getResult("IntDen",0);
		run("Clear Results");

		setResult("BranchNr",0,branch_nr);
		setResult("JunctionNr",0,junct_nr);
		setResult("EndpointNr",0,end_nr);
		setResult("TotalLength",0,branch_tot_length);
		setResult("AverageLength",0,branch_av_length);
		updateResults;	
			
		selectImage(idTagDup); close;
		selectImage(idTag); close;
		selectImage(idSkel); close;
		
		saveAs("Measurements", skeleton_results);
		erase(0);
		readout = 1;
	}	
	
	//	analyze nuclear rois
	if(segment_nuclei && File.exists(nuclei_roi_set))
	{
		run("Set Measurements...", "area mean standard modal min centroid center perimeter shape feret's integrated median skewness kurtosis redirect=None decimal=4");
		roiManager("Open",nuclei_roi_set);
		rmc = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
	
		sortResults(); //organise results per channel
				
		run("Select None");
		updateResults;
		saveAs("Measurements",nuclei_results);
		erase(0);
		readout = 1;
	}	
	
	
	if(!isOpen(mask)){
		newImage("Mask", "32-bit Black",image_width, image_height, 1); //	reference image for spot assignments 
		mask = getImageID;
	}
	
	//	rudimentary colocalization analysis by binary overlap of spot ROIs requires bianry masks of spots
	if(segment_spots && colocalize_spots && File.exists(spots_a_roi_set) && File.exists(spots_b_roi_set))
	{
		roiManager("reset");
		roiManager("Open",spots_a_roi_set);
		selectImage(mask);
		setForegroundColor(255,255,255);
		setSlice(1);
		roiManager("Fill");
		roiManager("reset");
		roiManager("Open",spots_b_roi_set);
		selectImage(mask);
		run("Add Slice");
		setSlice(2);
		roiManager("Fill");
		run("Select None");
		roiManager("reset");
	}
	
	//	analyze spot rois
	if(segment_spots && File.exists(spots_a_roi_set))
	{
		selectImage(mask); 
		ms = nSlices;
		run("Set Measurements...", "  area mean standard median min redirect=None decimal=4");
		roiManager("Open",spots_a_roi_set);
		spot_nr = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		
		// determine the colocalizing spots (one pixel overlap is sufficient)
		if(colocalize_spots && ms==2)
		{
			selectImage(mask); 
			setSlice(2);
			roiManager("Measure");
			overlaps = newArray(spot_nr);
			for(j=0;j<spot_nr;j++)
			{
				max = getResult("Max",j);
				if(max>0){overlaps[j]=1;}
			}	
			selectWindow("Results"); run("Close");
		}
		IJ.renameResults("Temp","Results");
		for(j=0;j<spot_nr;j++)
		{
			if(colocalize_spots && ms==2)setResult("Coloc",j,overlaps[j]);
		}
		updateResults;
		
		saveAs("Measurements",spots_a_results);
		erase(0);
		readout = 1;
	}
	if(segment_spots && File.exists(spots_b_roi_set))
	{
		selectImage(mask); 
		ms = nSlices;
		run("Set Measurements...", "  area mean standard median min redirect=None decimal=4");
		roiManager("Open",spots_b_roi_set);
		spot_nr = roiManager("count");
		selectImage(id);
		for(c=1;c<=channels;c++)
		{
			setSlice(c);
			roiManager("deselect");
			roiManager("Measure");
		}
		sortResults();
		IJ.renameResults("Results","Temp");
		
		// determine the colocalizing spots (one pixel overlap is sufficient)
		if(colocalize_spots && ms==2)
		{
			selectImage(mask); 
			setSlice(1);
			roiManager("Measure");
			overlaps = newArray(spot_nr);
			for(j=0;j<spot_nr;j++)
			{
				max = getResult("Max",j);
				if(max>0){overlaps[j]=1;}
			}	
			selectWindow("Results"); run("Close");
		}
		IJ.renameResults("Temp","Results");
		for(j=0;j<spot_nr;j++)
		{
			if(colocalize_spots && ms==2)setResult("Coloc",j,overlaps[j]);
		}
		updateResults;
		saveAs("Measurements",spots_b_results);
		erase(0);
		readout = 1;
	}
	if(isOpen(mask)){selectImage(mask); close;}
	return readout;
}

function summarizeResults(id)
{
	// 	open neurites results
	if(File.exists(search_results))
	{	
		run("Results... ", "open=["+search_results+"]");
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		for(r=0;r<1;r++)
		{
			for(s=0;s<resultLabels.length;s++)
			{
				selectImage(matrix);
				p = getPixel(s,r);
				setResult("Search_SC"+neurites_channel+"_"+resultLabels[s],r,p); // Label all cytoplasmic measured parameters with a "Cell" prefix
			}
		}
		updateResults;
		selectImage(matrix); close;
	}
	
	if(segment_neurites && File.exists(neurites_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+neurites_results+"]");
		resultLabels 	= getResultLabels();
		matrix 			= results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(r=0;r<1;r++)
		{
			for(s=0;s<resultLabels.length;s++)
			{
				selectImage(matrix);
				p = getPixel(s,r);
				setResult("Neurites_SC"+neurites_channel+"_"+resultLabels[s],r,p);
			}
		}
		updateResults;
		selectImage(matrix); close;
	}
	
	if(segment_neurites && File.exists(skeleton_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+skeleton_results+"]");
		resultLabels 	= getResultLabels();
		matrix 			= results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(r=0;r<1;r++)
		{
			for(s=0;s<resultLabels.length;s++)
			{
				selectImage(matrix);
				p = getPixel(s,r);
				setResult("Skeleton_SC"+neurites_channel+"_"+resultLabels[s],r,p);
			}
		}
		updateResults;
		selectImage(matrix); close;
	}
	
	//Summarise Nuclei
	if(File.exists(nuclei_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+nuclei_results+"]");
		// if distance measurement was done for spots add them here
		updateResults;
		nnr 			= nResults;
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		
		for(s=0;s<resultLabels.length;s++)
		{
			if(resultLabels[s] != "Coloc")
			{
				value=0;
				for(r=0;r<nnr;r++)
				{
					selectImage(matrix);
					p = getPixel(s,r);
					value += p;  
				}
				setResult("Nuclei_SC"+nuclei_channel+"_Nr",0,nnr);
				setResult("Nuclei_SC"+nuclei_channel+"_"+resultLabels[s]+"_Sum",0,value);              
				setResult("Nuclei_SC"+nuclei_channel+"_"+resultLabels[s]+"_Mean",0,value/nnr);
			}
		}
		selectImage(matrix); close;
		updateResults();
	}else if (segment_nuclei)
	{
		setResult("Nuclei_SC"+nuclei_channel+"_Nr",0,0);
	}
	
	//Summarise spots
	if(File.exists(spots_a_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+spots_a_results+"]");
		// if distance measurement was done for spots add them here
		updateResults;
		snr 			= nResults;
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(s=0;s<resultLabels.length;s++)
		{
			value=0;
			coloc=0;
			for(r=0;r<snr;r++)
			{
				selectImage(matrix);
				p = getPixel(s,r);
				if(resultLabels[s] != "Coloc")
				{	
					value += p;  
				}else if(p==1 && colocalize_spots){
					coloc=coloc+1;
				}
			}
			if(resultLabels[s] != "Coloc")
			{
				setResult("Spot_SC"+spots_a_channel+"_Nr",0,snr);
				setResult("Spot_SC"+spots_a_channel+"_"+resultLabels[s]+"_Sum",0,value);              
				setResult("Spot_SC"+spots_a_channel+"_"+resultLabels[s]+"_Mean",0,value/snr);
			}
			else{
				setResult("Spot_SC"+spots_a_channel+"_ColocNr",0,coloc);
				setResult("Spot_SC"+spots_a_channel+"_ColocPerc",0,coloc/snr*100);
			}
			
		}
		selectImage(matrix); close;
		updateResults();
	}
	
	if(File.exists(spots_b_results))
	{
		IJ.renameResults("Results","Temp");
		run("Results... ", "open=["+spots_b_results+"]");
		// if distance measurement was done for spots add them here
		updateResults;
		snr 			= nResults;
		resultLabels = getResultLabels();
		matrix = results2matrix(resultLabels);
		selectWindow("Results"); run("Close");
		IJ.renameResults("Temp","Results");
		for(s=0;s<resultLabels.length;s++)
		{
			
			value=0;
			coloc=0;
			for(r=0;r<snr;r++)
			{
				selectImage(matrix);
				p = getPixel(s,r);
				if(resultLabels[s] != "Coloc")
				{	
					value += p;  
				}else if(p==1){
					coloc=coloc+1;
				}
			}
			if(resultLabels[s] != "Coloc")
			{
				setResult("Spot_SC"+spots_b_channel+"_Nr",0,snr);
				setResult("Spot_SC"+spots_b_channel+"_"+resultLabels[s]+"_Sum",0,value);              
				setResult("Spot_SC"+spots_b_channel+"_"+resultLabels[s]+"_Mean",0,value/snr);
			}
			else{
				setResult("Spot_SC"+spots_b_channel+"_ColocNr",0,coloc);
				setResult("Spot_SC"+spots_b_channel+"_ColocPerc",0,coloc/snr*100);
			}
			
		}
		selectImage(matrix); close;
		updateResults();
	}
	
	selectImage(id);
	if(segment_spots && colocalize_spots & File.exists(search_roi_set)){
		roiManager("reset");
		roiManager("open", search_roi_set);
		roiManager("select", 0);
		pearson=pearsonCorrelation(id,spots_a_channel,spots_b_channel);
		setResult("Spots_SC"+spots_a_channel+"_SC"+spots_b_channel+"_PearsonCorrelation", 0, pearson);
		updateResults();
	}
	roiManager("reset");
	selectWindow("Results"); saveAs("Measurements",results);
}

function sortResults()
{
	resultLabels = getResultLabels();
	matrix = results2matrix(resultLabels);
	matrix2results(matrix,resultLabels,channels);
}

function getResultLabels()
{
	selectWindow("Results");
	ls 				= split(getInfo(),'\n');
	rr 				= split(ls[0],'\t'); 
	nparams 		= rr.length-1;			
	resultLabels 	= newArray(nparams);
	for(j=1;j<=nparams;j++){resultLabels[j-1]=rr[j];}
	return resultLabels;
}

function results2matrix(resultLabels)
{
	h = nResults;
	w = resultLabels.length;
	newImage("Matrix", "32-bit Black",w, h, 1);
	matrix = getImageID;
	for(j=0;j<w;j++)
	{
		for(r=0;r<h;r++)
		{
			v = getResult(resultLabels[j],r);
			selectImage(matrix);
			setPixel(j,r,v);
		}
	}
	run("Clear Results");
	return matrix;
}

function matrix2results(matrix,resultLabels,channels)
{
	selectImage(matrix);
	w = getWidth;
	h = getHeight;
	for(c=0;c<channels;c++)
	{
		start = c*h/channels;
		end = c*h/channels+h/channels;
		for(k=0;k<w;k++)
		{
			for(j=start;j<end;j++)
			{
				selectImage(matrix);
				p = getPixel(k,j);
				setResult(resultLabels[k]+"_MC"+c+1,j-start,p); // MC for measurement channel
			}
		}
	}
	selectImage(matrix); close;
	updateResults;
}

function matrix2resultsConcatenate(matrix,resultLabels)
{
	selectImage(matrix);
	w = getWidth;
	h = getHeight;
	start = nResults;
	for(k = 0;k < w; k++)
	{
		for(j = 0; j < h; j++)
		{
			selectImage(matrix);
			p = getPixel(k,j);
			setResult(resultLabels[k],j+start,p); 
		}
	}
	selectImage(matrix); close;
	updateResults;
}

function toggleOverlay()
{	
	run("Select None"); 
	roiManager("deselect");
	roiManager("Show All without labels");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}
	
function createOverlay(names)
{
	setForegroundColor(25, 25, 25);
	fields = names.length;
	print(fields,"images");
	for(i=0;i<fields;i++)
	{
		prefix = names[i];
		file = prefix+suffix;
		setFileNames(prefix);
		print(i+1,"/",fields,":",prefix);
		path = dir+file;
		run("Bio-Formats Importer", "open=["+path+"] color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		//open(path);
		id = getImageID;
		Stack.getDimensions(w,h,channels,slices,frames); 
		if(!Stack.isHyperStack && channels == 1)
		{
			channels = slices;
			run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices=1 frames=1 display=Composite");
		}
		id = getImageID;

		if(segment_nuclei){
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(nuclei_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",nuclei_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
		
		if(segment_neurites){
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(neurites_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",neurites_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
			
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(skeleton_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",skeleton_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
		
		selectImage(id);
		setSlice(nSlices);
		run("Add Slice","add=channel");
		if(File.exists(search_roi_set))
		{	
			selectImage(id);
			setSlice(nSlices);
			roiManager("Open",search_roi_set);
			roiManager("deselect");
			roiManager("Fill");
			roiManager("reset");
		}
		
		if(segment_spots && spots_a_channel>0)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_a_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_a_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
		
		if(segment_spots && spots_b_channel>0)
		{
			selectImage(id);
			setSlice(nSlices);
			run("Add Slice","add=channel");
			if(File.exists(spots_b_roi_set))
			{	
				selectImage(id);
				setSlice(nSlices);
				roiManager("Open",spots_b_roi_set);
				roiManager("deselect");
				roiManager("Fill");
				roiManager("reset");
			}
		}
	}
	run("Concatenate...", "all_open title=[Concatenated Stacks]");
	Stack.getDimensions(w,h,newchannels,slices,frames);
	for(c=1;c<=channels;c++){Stack.setChannel(c);Stack.setFrame(round(frames/2));resetMinAndMax;}
	range = pow(2,bitDepth);
	for(c=channels+1;c<=newchannels;c++){Stack.setChannel(c);setMinAndMax(0,range/2);}
	
	idConc1=getImageID();
	if(segment_spots){
		run("Duplicate...", "duplicate");
		idConc2=getImageID();
	}
	
	selectImage(idConc1);
	Stack.getDimensions(w,h,newchannels,slices,frames);
	if(segment_spots && spots_b_channel >0){
		Stack.setChannel(newchannels);
		run("Delete Slice", "delete=channel");
	}
	Stack.getDimensions(w,h,newchannels,slices,frames);
	if(segment_spots && spots_a_channel >0){
		Stack.setChannel(newchannels);
		run("Delete Slice", "delete=channel");
	}	
	run("Make Composite");
	
	if(segment_spots){
		selectImage(idConc2);
		if(segment_nuclei){
			Stack.setChannel(channels+1);
			run("Delete Slice", "delete=channel");
		}
		if(segment_neurites){
			Stack.setChannel(channels+1);
			run("Delete Slice", "delete=channel");
			Stack.setChannel(channels+1);
			run("Delete Slice", "delete=channel");
		}	
		run("Make Composite");		
	}
	
}

function concatenateResults()
{
	prefixes  = scanFiles();
	index = 0;
	names = newArray(0);
	for(i = 0; i < prefixes.length; i++)
	{
		print(i+1,"/",prefixes.length,"...");
		results = output_dir+prefixes[i]+"_Summary.txt";
		if(File.exists(results)){
			run("Results... ","open=["+results+"]");
			nr = nResults; print("...",nr,"results");
			for(r = index; r < index+nr; r++)names = Array.concat(names,prefixes[i]); 
			if(i > 0) 
			{	
				labels = getResultLabels();
				matrix = results2matrix(labels);
				selectWindow("Results"); run("Close");
				IJ.renameResults("Summary","Results");
				matrix2resultsConcatenate(matrix,labels);
			}
			selectWindow("Results"); 
			IJ.renameResults("Summary");
		}
	}
	IJ.renameResults("Summary","Results");
	for(i = 0; i < names.length; i++)setResult("Label",i,names[i]);
	updateResults;
}

function pearsonCorrelation (id,ch1,ch2){
	selectImage(id);
	Stack.setChannel(ch1);run("Duplicate...","title=C1");C1id=getImageID;
	selectImage(C1id);if(bitDepth<32)run("32-bit");run("Select None");
	selectImage(id);
	Stack.setChannel(ch2);run("Duplicate...","title=C2");C2id=getImageID;
	selectImage(C2id);if(bitDepth<32)run("32-bit");run("Select None");
	
	selectImage(C1id);run("Duplicate...","title=C1n");
	roiManager("select", 0);
	getRawStatistics(n,m1,min1,max1);
	run("Subtract...", "value="+m1);
	
	
	run("Duplicate...","title=C1ns");
	run("Square");
	roiManager("select", 0);
	getRawStatistics(n,m1ns);
	sns1=n*m1ns;		
	selectImage(C2id);run("Duplicate...","title=C2n");
	roiManager("select", 0);
	getRawStatistics(n,m2,min2,max2);
	run("Subtract...", "value="+m2);
	
	run("Duplicate...","title=C2ns");
	run("Square");
	roiManager("select", 0);
	getRawStatistics(n,m2ns);sns2=n*m2ns;
	
	
	run("Image Calculator...", "image1=C1n operation=Multiply image2=C2n create 32-bit"); 		
	selectWindow("Result of C1n");
	roiManager("select", 0);
	getRawStatistics(n,m12n);
	s12n=n*m12n;
	rename("PDM");
	selectWindow("C1ns");close();
	selectWindow("C2ns");close();
	selectWindow("C1n");close();
	selectWindow("C2n");close();
	selectWindow("PDM");close();
	selectImage(C1id);close();
	selectImage(C2id);close();
	
	pearson=s12n/sqrt(sns1*sns2);
	
	return pearson;	
}
