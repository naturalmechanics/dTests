module analysisEngine;

import field;
import neuron;

import std.stdio;
import std.conv;
import std.format;
import std.math;
import std.algorithm;
import std.net.curl;
import std.json;
import dlib.image;
import std.path;
import std.array;
import std.net.curl;
import core.stdc.stdlib;
import std.datetime;
import std.file;


class geoEngine  {

	this(){}
	~this(){}


public :

	void * dataSet ;
	int[] convexHull;	
	string fileName;
	
	
	double maxLineDist			= 30;
	double minSectorLength 		= 10;
	double maxDeviation_ofAngle = 15.00 * PI / 180.00;
	double maxLineLengthOvershotRatio = 0.1;
	double maxLineOutlierCount  = 5;
	double minLineLength		= 8;
	int    minLineCount         = 10;
	int    minLineCount_forTrim = 2;
	int minLineCount_nearVertex = 5;
	double crad                 = 40;
	double crad2				= 40;
	double minCrossLineLength   = 10;
	double maxTurnRadius		= 15;
	ulong maxMeridianSections   = 50;
	
	double angularResolution    = 0.1;
	double minSpatialResolution = 0.2;
	double lineLengthScoreAmplifier		= 4;
	double lineLengthScoreWeight        = 10;
	double lineLengthScoreBias			= 0.2;
	double stripeThinnessScoreAmplifier = 100;
	double stripeThinnessScoreWeight	= 12;
	double stripeThinnessScoreBias 		= 2;
	
	
	string default_maxDrivingFluctuation = "10";
	int maxMergeAttempts = 150;
	int maxBoundaryReductionSteps = 50000;
	
	double transportCircleRadius= 20;
	double minPointDistance = 2.0;
	double max_interTankdistance = 40;
	double transport_fieldThreshhold = 250;
	int    transport_xgressLength= 15;
	double loadingPointThreshold = 10.0;
	int holeClosingThreshhold 	 = 2;
	
	double workWidth            = 580;
	double geoFenceRadius       = 200;
	double[] geoFenceCoordinates= new double [] (0);
	double geoFenceRadiusThreshhold = 50;
	
	int    trendAveragingLength =4;
	double parallelThreshholdResolution = 0.1;
	double maximumAllowedNonParallelism = 80;
	double minOverlap           = 0.5;
	
	double scoreCutOff          = 1.7;
	
	enum possibleDatum : int  {WGS1984};
	
	possibleDatum actualDatum = possibleDatum.WGS1984;
	
	
	
	void set_dataType(int typ) {															// JUST a set function. no corresponding get function
		switch (typ)  {
		
			case 1:
			case 2:
				this.dataType_current = 2;
			break;
			default :
			break;

		}
	}
	
	void calculate_drawingParams() {														// no arguments, 
																				// this will automatically take the main dataset

		switch (this.dataType_current) {

			case 1 :															// two arrays, with X and Y
			case 2 :															// geodata 
				
				double latMax  = -90;
				double latMin  =  90;
				double lonMax  =-180;
				double lonMin  = 180;
				
				int leftMostPos;
				
				field.rawData [] * rd ;
				rd = cast (field.rawData [] *)  dataSet;						// rawData rd is the complete dataset of all points

				foreach (point; *rd) {											// each point has 
																				// one lat, one lon, one time string, one satcount.
																				// point type is field.rawData
					if (point.lat > latMax)		latMax = point.lat;
					if (point.lat < latMin)		latMin = point.lat;
					if (point.lon > lonMax)		lonMax = point.lon;
					if (point.lon < lonMin)		lonMin = point.lon;
																				// writeln( "Boundary found : " ~ format("%.*g", 18, latMax) ~ "; " ~ format("%.*g", 18, lonMax) ~ "; " ~ format("%.*g", 18,latMin) ~ "; " ~  format("%.*g", 18,lonMin) ~ "; ");
				}
				
				double latDiff = latMax - latMin;
				double lonDiff = lonMax - lonMin;
				
				this.latMin_global = latMin;
				this.latMax_global = latMax;
				this.lonMin_global = lonMin;
				this.lonMax_global = lonMax;
				
				this.latdiff_global = latDiff;
				this.londiff_global = lonDiff;
				
			break;
			default:
			break;
		}
	}
	
	void create_map_fromRawData_inDLIB (double width, double height) {
	
		if(width == -1)  {width  = this.londiff_global*1.1;}
		if(height == -1) {height = this.latdiff_global*1.1;}								// automatic mode engaged
																							writeln([width,height]);
		double multiplier = 1;
		
		
		if (width < 0.1 || height < 0.1 ) {													writeln("width and height too small need rescale..");
			
			multiplier = 800000;															write("image width is : ");writeln(width);
			width = width * multiplier *cos(this.latMin_global * 3.141592 / 180.00);		write("image height is: ");writeln(height);
			height = height * multiplier;
			
		}
		else if (width < 1 || height < 1 ) {
			
			multiplier = 800000;
			width = width * multiplier *cos(this.latMin_global * 3.141592 / 180.00);
			height = height * multiplier;
			
		}
		
		
		
		this.multiplier_global = multiplier;
		
																							write("diffs are: ");writeln([this.londiff_global, this.latdiff_global]);
																							write("recheck dimensions - those are: ");writeln([width,height]);
		auto image = new Image!(PixelFormat.RGB8)(to!int(width), to!int(height));	
																							write("after conversion to int : ");writeln([to!int(width), to!int(height)]);
// 		foreach(y; 0 .. image.height){
// 			foreach(x; 0 .. image.width){
// 				image[x,y] = Color4f(0.1,0.1,0.1,1);
// 			}
// 		}
		
		double offsetY = londiff_global * 0.05 * multiplier;
		double offsetX = latdiff_global * 0.05 * multiplier;
		
		field.rawData [] * rd ;
		rd = cast (field.rawData [] *)  dataSet;											// rawData rd is the complete dataset of all points

		foreach (point; *rd) {																// each point has 
																							// one lat, one lon, one time string, one satcount.
																							// and one UNIQUE ID
																							// point type is field.rawData
			int mapy = to!int(multiplier * (this.latMax_global - point.lat) + offsetX);
			int mapx = to!int((multiplier * (point.lon - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
			
			image[mapx,mapy] = Color4f(1,1,1,1);											// write("drawing");writeln([mapx,mapy]);
			
		}
				
				
		this.map = image;
		this.mapCopy = image;
		
	}

	void analyze_cultivation() {
		
// 		get_allLines();																	writeln("all lines found : Check color Cyan");
// 		guess_fieldLines();																writeln("all lines found : Check color Green");
// 		combine_fieldLines();
// 		guess_boundaries();																writeln("boundaries found: Check color Yellow");
// 		trim_fieldLines();
// 		
// 		get_rawBoundaries();
// 		
// 		get_bridgeLines();
// 		
// 		merge_fields();
// 		get_redLines();
// 		get_geoFenceEvents();
// 		trim_fieldLines();																writeln("trimming found  : Check color Red");
// 		submit_results();
 		
		set_lineDetectorNeurons();														// setting neurons 
		/+fire_lineDetectorNeurons();
		findSpikes_from_lineDetectorNeurons();+/
	}
	
	void analyze_tillage() {
	
// 		double offsetY = londiff_global * 0.05 * this.multiplier_global;
// 		double offsetX = latdiff_global * 0.05 * this.multiplier_global;
// 
// 		field.rawData [] * rd ;
// 		rd = cast (field.rawData [] *)  dataSet;
// 
// 		get_allLines();
// 		guess_fieldLines();
// 		guess_boundaries();
// 		
// 		trim_fieldLines();
// 		
// 		get_rawBoundaries();
// 		
// 		get_bridgeLines();
// 		
// 		merge_fields();
// 		
// 		
// 		
// 		
// 		for(int i = 0; i < this.fieldBoundariesRaw.length; i++) {
// 		
// 			auto F = this.fieldBoundariesRaw[i];								writeln(F);
// 			
// 			for (int j = 0; j < F.length-1; j ++){
// 			
// 				auto v0 = F[j];
// 				auto v1 = F[j+1];
// 				
// 				double y1 = to!int(this.multiplier_global * (this.latMax_global - v0[1]) + offsetX);
// 				double y2 = to!int(this.multiplier_global * (this.latMax_global - v1[1]) + offsetX);
// 				double x1 = to!int((this.multiplier_global * (v0[0] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
// 				double x2 = to!int((this.multiplier_global * (v1[0] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
// 				auto RED = Color4f(250.00/255.00,50.00/255.00,50.00/255.00);							//writeln(0);
// 				drawLine_png(this.map,RED,to!int(x1),to!int(y1),to!int(x2),to!int(y2));
// 				
// 			
// 			}
// 			
// 		}
// 		
// 		for( int i = 0; i < this.baseFields.length; i++) {
// 		
// 			auto oneField = this.baseFields[i];
// 			double tl = 0;
// 			int gstep;
// 			for (int j = 0; j < oneField.length; j++) {
// 
// 				auto oneLine = oneField[j];
// 				
// 				for(gstep = this.allLines[oneLine][0]; gstep < this.allLines[oneLine][1]; gstep++) {
// 	
// 					auto P1 = (*rd)[gstep];
// 					auto P2 = (*rd)[gstep +1];
// 
// 					double y1 = to!int(this.multiplier_global * (this.latMax_global - P1.lat) + offsetX);
// 					double y2 = to!int(this.multiplier_global * (this.latMax_global - P2.lat) + offsetX);
// 					double x1 = to!int((this.multiplier_global * (P1.lon - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
// 					double x2 = to!int((this.multiplier_global * (P2.lon - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
// 					auto GREEN = Color4f(20.00/255.00,250.00/255.00,250.00/255.00);							//writeln(0);
// 					drawLine_png(this.map,GREEN,to!int(x1),to!int(y1),to!int(x2),to!int(y2));
// 																	// array for one line close
// 				}
// 				
// 				
// 				
// 				
// 			}
// 		}
// 		
// 		//guess_crossLines();
// 		merge_fields();
// 		
// 		draw_map(); 
// 		// readln();
// 		get_redLines();
// 		get_geoFenceEvents();
//  		submit_results();

//		set_lineDetectorNeurons();															writeln("all neurones set");
//		fire_lineDetectorNeurons();															writeln("all neurones fired");
//		findSpikes_from_lineDetectorNeurons();												writeln("SCORES found");
		//get_allLines();
	}
	
	
	
	
	
	
	
	
	void set_lineDetectorNeurons() {
	
		double left    = this.lonMin_global;
		double right   = this.lonMax_global;
		double Top     = this.latMax_global;
		double Bottom  = this.lonMax_global;
		
		
		string name = baseName(stripExtension(this.fileName));								// writeln("file name found" ~ name);
		auto parts  = split(name, "-");														// found the devid and pathname
		auto devID  = parts[0];																// we will need the dev ID to search memory
		
		double offsetY = londiff_global * 0.05 * this.multiplier_global;
		double offsetX = latdiff_global * 0.05 * this.multiplier_global;
		
		for (double angle = 0; angle < PI ; angle = angle + this.angularResolution) {		//write("Testing angle : ");writeln(angle);
																			
			// double extent = 0; // get_fieldExtent(left, right, top, bottom, angle) ;
			
			double max_neuronAreaWidth = to!double(get_fromMemory(devID, "maxDrivingFluctuation"));
			
			for (double width = this.minSpatialResolution ; width <= max_neuronAreaWidth ; width = width + this.minSpatialResolution) {
				
				if (angle == 0 ) {
					
					neuron.edgeDetectorNeuron eDN ;
				
					eDN.angle = angle;
					eDN.width = width;
					
					auto lastP1 = [this.latMax_global, this.lonMax_global];					// start at top right
					auto lastP2 = [this.latMax_global, this.lonMin_global];
					
					double iw3 = width;
					double totLat = this.latMax_global - this.latMin_global;
					double latD   = calculate_geoDistance (this.latMax_global, this.latMin_global, this.lonMin_global, this.lonMin_global); // lonmin because start at left
					int latsegs	  = to!int(floor(latD / iw3));
					
					double latStep= totLat / latsegs; 
					
					bool breakFlag = false;
					
					while( ! breakFlag) {
					
						auto tent_nextP3 = [lastP2[0] -latStep, lastP2[1]];
						auto tent_nextP4 = [lastP1[0] -latStep, lastP1[1]];
						
						if (tent_nextP3[0] < this.latMin_global &&  tent_nextP4[0] < this.latMin_global) {
						
							breakFlag = true;
							
							tent_nextP3 = [lastP2[0], lastP2[1]];
							tent_nextP4 = [lastP2[0], lastP2[1]];
							
						}
						
						auto nextP3 = tent_nextP3;
						auto nextP4 = tent_nextP4;
						
						neuron.edgeDetectorNeuron_area na;
						na.v1 = lastP1;
						na.v2 = lastP2;
						
						na.v3 = nextP3;
						na.v4 = nextP4;
						
						lastP1 = nextP4;													// the right hand side of the polygon will become the left hand side of the next polygon
						lastP2 = nextP3;													// this is why we flip the sequence : 3 ----> 2 , 4 ----> 1 
						
						eDN.stripe = na;	
						
						this.edgeDetectors ~= eDN;
						
						
						
					
					}
				} else if (angle > 0 && angle < PI / 2.00) {
																							// writeln ( "get leftmost intersectpoints");
																							// writeln ( "in this case starts at the top left");
																							// writeln ( "vertices are counted frop top left, counterclokwise");
					
					neuron.edgeDetectorNeuron eDN ;
				
					eDN.angle = angle;
					eDN.width = width;
					
					
					auto lastP1 = [this.latMax_global, this.lonMin_global];
					auto lastP2 = [this.latMax_global, this.lonMin_global];
					
					double iw4 = width / sin( angle);
					double iw3 = width / cos (angle);										//writeln(iw3);writeln(width);
						
					double totLon = this.lonMax_global - this.lonMin_global;
					double lonD   = calculate_geoDistance (this.latMax_global, this.latMax_global, this.lonMin_global, this.lonMax_global);	// latmax, becasue start at top
					int lonsegs	  = to!int(floor(lonD / iw4));
					
					if (lonsegs == 0) continue;
					
					double lonStep= totLon / lonsegs; 
					
						
					double totLat = this.latMax_global - this.latMin_global;
					double latD   = calculate_geoDistance (this.latMax_global, this.latMin_global, this.lonMin_global, this.lonMin_global); // lonmin because start at left
					int latsegs	  = to!int(floor(latD / iw3));								 //writeln(latsegs);
					
					if (latsegs == 0) continue;
					
					double latStep= totLat / latsegs; 
					
					bool breakFlag = false;
					
					while( ! breakFlag) {													// writeln("starting ...");
						
						auto tent_nextP3 = [lastP2[0] -latStep, lastP2[1]];
						auto tent_nextP4 = [lastP1[0], lastP1[1] + lonStep];					
																							// writeln(toStringLikeInCSharp(tent_nextP3[0]));
																							// writeln(toStringLikeInCSharp(tent_nextP3[1]));
																							// writeln(toStringLikeInCSharp(tent_nextP4[0]));
																							// writeln(toStringLikeInCSharp(tent_nextP4[1]));
						
						if (tent_nextP3[1] >= this.lonMax_global &&  tent_nextP4[0] <= this.latMin_global) { // last point, will break after this
																							//writeln("ready to break");
							breakFlag = true;
							tent_nextP3 = [this.latMin_global, this.lonMax_global];
							tent_nextP4 = [this.latMin_global, this.lonMax_global];
						} else {
						
							if (tent_nextP3 [0] <= this.latMin_global) { 					// crossed the bottom left limit, will increase lon instead
																							// writeln("bottom reached");
								tent_nextP3 = [ this.latMin_global , min(lastP2[1] + lonStep, this.lonMax_global)];
							} 
							
							if (tent_nextP4 [1] >= this.lonMax_global) {					// crossed the top right boundary
																							// writeln("right end reached");
								tent_nextP4 = [max(lastP1[0] - latStep, this.latMin_global),  this.lonMax_global ];
								
							}
						}
						
																							// writeln("while continues");
																							// writeln(toStringLikeInCSharp(tent_nextP3[0]));
																							// write("--->");write(toStringLikeInCSharp(tent_nextP3[1]));
																							// write(";");writeln(toStringLikeInCSharp(this.lonMax_global));
																							// write("--->");write(toStringLikeInCSharp(tent_nextP4[0]));
																							// write(";");writeln(toStringLikeInCSharp(this.latMin_global));
																							// writeln(toStringLikeInCSharp(tent_nextP4[1]));
																				
																				
						auto nextP3 = tent_nextP3;
						auto nextP4 = tent_nextP4;
						
						neuron.edgeDetectorNeuron_area na;
						na.v1 = lastP1;
						na.v2 = lastP2;
						
						na.v3 = nextP3;
						na.v4 = nextP4;
						
						lastP1 = nextP4;													// the right hand side of the polygon will become the left hand side of the next polygon
						lastP2 = nextP3;													// this is why we flip the sequence : 3 ----> 2 , 4 ----> 1 
						
						eDN.stripe = na;													// writeln(eDN);
						
						this.edgeDetectors ~= eDN;
						
							
					}
					
					
					
					
				} else if (angle == PI / 2.00) {
				
					neuron.edgeDetectorNeuron eDN ;
				
					eDN.angle = angle;
					eDN.width = width;
				
					auto lastP1 = [this.latMax_global, this.lonMin_global];					// start at top Left
					auto lastP2 = [this.latMin_global, this.lonMin_global];
					
					double iw4 = width;
					
					double totLon = this.lonMax_global - this.lonMin_global;
					double lonD   = calculate_geoDistance (this.latMax_global, this.latMax_global, this.lonMin_global, this.lonMax_global);	// latmax, becasue start at top
					int lonsegs	  = to!int(floor(lonD / iw4));
					
					double lonStep= totLon / lonsegs; 
					
					bool breakFlag = false;
					
					while( ! breakFlag) {
					
						auto tent_nextP3 = [lastP2[0] , lastP2[1] + lonStep];
						auto tent_nextP4 = [lastP1[0] , lastP1[1] + lonStep];
						
						if (tent_nextP3[1] > this.lonMax_global &&  tent_nextP4[1] > this.lonMax_global) {
						
							breakFlag = true;
							
							tent_nextP3 = [lastP2[0], lastP2[1]];
							tent_nextP4 = [lastP2[0], lastP2[1]];
							
						}
						
						auto nextP3 = tent_nextP3;
						auto nextP4 = tent_nextP4;
						
						neuron.edgeDetectorNeuron_area na;
						na.v1 = lastP1;
						na.v2 = lastP2;
						
						na.v3 = nextP3;
						na.v4 = nextP4;
						
						lastP1 = nextP4;													// the right hand side of the polygon will become the left hand side of the next polygon
						lastP2 = nextP3;													// this is why we flip the sequence : 3 ----> 2 , 4 ----> 1 
						
						eDN.stripe = na;	
						
						this.edgeDetectors ~= eDN;
						
					
					}
					
				
				} else if (angle > PI / 2.00) {
																							// writeln ( "get leftmost intersectpoints");
																							// writeln ( "in this case starts at the bottom left");
																							// writeln ( "vertices are counted frop top left, counterclokwise");
					
					
					neuron.edgeDetectorNeuron eDN ;
				
					eDN.angle = angle;
					eDN.width = width;
					
					
					auto lastP1 = [this.latMin_global, this.lonMin_global];
					auto lastP2 = [this.latMin_global, this.lonMin_global];
					
					double iw4 = width / sin (angle - PI / 2.00);
					double iw3 = width / cos (angle - PI / 2.00);
						
					double totLon = this.lonMax_global - this.lonMin_global;
					double lonD   = calculate_geoDistance (this.latMax_global, this.latMax_global, this.lonMin_global, this.lonMax_global);	// latmax, becasue start at top
					int lonsegs	  = to!int(floor(lonD / iw3));
					
					if (lonsegs == 0) continue;
					
					double lonStep= totLon / lonsegs; 
					
						
					double totLat = this.latMax_global - this.latMin_global;
					double latD   = calculate_geoDistance (this.latMax_global, this.latMin_global, this.lonMin_global, this.lonMin_global); // lonmin because start at left
					int latsegs	  = to!int(floor(latD / iw4));
					
					if (lonsegs == 0) continue;
					
					double latStep= totLat / latsegs; 
					
					bool breakFlag = false;
					
					while( ! breakFlag) {
						
						auto tent_nextP3 = [lastP2[0], lastP2[1] + lonStep];
						auto tent_nextP4 = [lastP1[0] + latStep, lastP1[1]];
						
						if (tent_nextP3[0] >= this.latMax_global &&  tent_nextP4[1] >= this.lonMax_global) { // last point, will break after this
							breakFlag = true;												//readln();
							tent_nextP3 = [this.latMax_global, this.lonMax_global];
							tent_nextP4 = [this.latMax_global, this.lonMax_global];
						} else {
						
							if (tent_nextP3 [1] >= this.lonMax_global) { 					// crossed the bottom left limit, will increase lon instead
								tent_nextP3 = [min(lastP2[0] + latStep, this.latMax_global) , this.lonMax_global];
							} 
							
							if (tent_nextP4 [0] >= this.latMax_global) {					// crossed the top right boundary
								tent_nextP4 = [this.latMax_global , min(lastP1[1] + lonStep, this.lonMax_global) ];
							}
						
						}
						
						auto nextP3 = tent_nextP3;
						auto nextP4 = tent_nextP4;
						
						neuron.edgeDetectorNeuron_area na;
						na.v1 = lastP1;
						na.v2 = lastP2;
						
						na.v3 = nextP3;
						na.v4 = nextP4;
						
						lastP1 = nextP4;													// the right hand side of the polygon will become the left hand side of the next polygon
						lastP2 = nextP3;													// this is why we flip the sequence : 3 ----> 2 , 4 ----> 1 
						
						eDN.stripe = na;										
						
						this.edgeDetectors ~= eDN;											// writeln("drawing");
					
					
						auto ANG_pi_over_2_to_pi = Color4f((50.00*angle)/255.00,(20.00*eDN.width)/255.00,20.00/255.00);	//writeln(0);
						
						auto yy11 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v1[0]) + offsetX);
						auto yy12 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v2[0]) + offsetX);
						auto xx11 = to!int((this.multiplier_global * (eDN.stripe.v1[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						auto xx12 = to!int((this.multiplier_global * (eDN.stripe.v2[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						drawLine_png(this.map,ANG_pi_over_2_to_pi,to!int(xx11),to!int(yy11),to!int(xx12),to!int(yy12));
						
						auto yy21 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v2[0]) + offsetX);
						auto yy22 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v3[0]) + offsetX);
						auto xx21 = to!int((this.multiplier_global * (eDN.stripe.v2[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						auto xx22 = to!int((this.multiplier_global * (eDN.stripe.v3[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						drawLine_png(this.map,ANG_pi_over_2_to_pi,to!int(xx21),to!int(yy21),to!int(xx22),to!int(yy22));
						
						auto yy31 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v3[0]) + offsetX); 
																								// writeln(this.multiplier_global * (this.latMax_global - eDN.stripe.v4[0]) + offsetX); 
																								// writeln( (this.latMax_global - eDN.stripe.v4[0]) + offsetX);
						auto yy32 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v4[0]) + offsetX);
						auto xx31 = to!int((this.multiplier_global * (eDN.stripe.v3[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						auto xx32 = to!int((this.multiplier_global * (eDN.stripe.v4[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						drawLine_png(this.map,ANG_pi_over_2_to_pi,to!int(xx31),to!int(yy31),to!int(xx32),to!int(yy32));
						
						auto yy41 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v4[0]) + offsetX);
						auto yy42 = to!int(this.multiplier_global * (this.latMax_global - eDN.stripe.v1[0]) + offsetX);
						auto xx41 = to!int((this.multiplier_global * (eDN.stripe.v4[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						auto xx42 = to!int((this.multiplier_global * (eDN.stripe.v1[1] - this.lonMin_global) + offsetY) * cos(this.latMin_global * 3.141592 / 180.00));
						drawLine_png(this.map,ANG_pi_over_2_to_pi,to!int(xx41),to!int(yy41),to!int(xx42),to!int(yy42));
							
						
					}																		// writeln([xx41, xx42, yy41, yy42]);
					
					
				}
				
			}
		
		}																					writeln(this.edgeDetectors.length); draw_map();
		
		
	}
	
	void set_lineDetectorNeurons(int resolution, int orientation) {
	
	}
	
	void set_lineDetectorNeurons(int resolution, int orientation, int nums) {
	
	}
	
	void draw_map() { 
																							writeln("waiting to finish drawing"); 
		savePNG(this.map, fileName~".png"); 												writeln("check map: "~fileName ~".png"); 
	}
	
	
	
	void drawLine_png(Image!(PixelFormat.RGB8) img, Color4f color, int x1, int y1, int x2, int y2) {
		int dx = x2 - x1;														
		int ix = (dx > 0) - (dx < 0);											
		int dx2 = std.math.abs(dx) * 2;											
		int dy = y2 - y1;														
		int iy = (dy > 0) - (dy < 0);											
		int dy2 = std.math.abs(dy) * 2;
		img[x1, y1] = color;
																				
		if (dx2 >= dy2)
		{																		
			int error = dy2 - (dx2 / 2);
			while (x1 != x2)
			{
				if (error >= 0 && (error || (ix > 0)))
				{
					error -= dx2;
					y1 += iy;
				}

				error += dy2;
				x1 += ix;
				img[x1, y1] = color;
			}
		} else {																		
			int error = dx2 - (dy2 / 2);
			while (y1 != y2)
			{
				if (error >= 0 && (error || (iy > 0)))
				{
					error -= dy2;
					x1 += ix;
				}

				error += dx2;
				y1 += iy;
				img[x1, y1] = color;
			}
		}
	}
	
private :

	short dataType_current;
	int drawH, drawW;
	int convexHullLength;
	int[] semiMajorAxis;
	int[] semiMinorAxis;
	
	double latdiff_global;
	double londiff_global;
	
	double latMax_global;
	double latMin_global;
	double lonMax_global;
	double lonMin_global;
	double multiplier_global;
	
	int transport_startEpicenter;
	int transport_endEpicenter;
	int xGressCount;
	
	neuron.edgeDetectorNeuron[] edgeDetectors ;
	
	
	int [][] baseFields = new int[][] (0,0);
	int [][] plainFields= new int[][] (0,0);
	int [][] fieldBoundaries = new int [][] (0,0);
	double [][][] fieldBoundariesRaw = new double [][][] (0,0,0);
	double [][] emergency = new double  [][] (0,0);
	int [] allRedLines = new int[] (0);
	
	int [] ingressPoints = new int[](0);
	int [] outgressPoints = new int[](0);
	
	
	double[][][] innerTriangles = new double [][][] (0,0,0);
	int [] candidatePoints = new int [] (0);
	int[][] bestFieldLines = new int [][] (0,0);
	int [] fieldIndices = new int [] (0);
	int [][] allLines = new int [][] (0,0);
	int [] goodlines = new int [] (0);
	int [][][] fieldCorrections = new int [][][] (0,0,0);
	
	int[][] knotPointIndices = new int [][] (0,0);
	int mostLikelyLocation_ofFarm = -1;
	int [] mergedLocation_ofFarm = new int[] (0);
	double [] farmCenter = new double [] (0);
	int firstLoadingPoint = -1;
	
	int [][] skiplines = new int[][] (0,0);
	int [] notRedLinesGlobal = new int [](0);
	
	auto map = new Image!(PixelFormat.RGB8)(1,1)  ;
	auto mapCopy =  new Image!(PixelFormat.RGB8)(1,1)  ;
	auto convexHullDiagonals_ofMap = new int[][] (0,0);
	
	
	double fieldSizeNormalizer = 8192;
	
	
	string jsonstr;
	
	int[][] estimatedFields;
	int [][] greenLimits = new int [][] (0,0);
	
	
	
	
	
	
	
	double calculate_geoDistance(double lat1, double lat2, double lon1, double lon2) {
	
		double a = pow(( sin ( ( lat1 * 3.141592/180.00 - lat2 * 3.141592 / 180.00 )  / 2.00) ) , 2);
		double b = cos( lat1 * 3.141592 / 180.00 ) * cos(lat2 * 3.141592 / 180.00);
		double c = pow(( sin ( ( lon1 * 3.141592 / 180.00 - lon2 * 3.141592 / 180.00) / 2.00)  ) ,2);
		
		double dist = 2 * 6400000 * asin (sqrt ( a+ b * c )	);
		return dist;

	}

	string get_fromMemory(string key, string source) {
	
		string res = this.default_maxDrivingFluctuation;
	
	return res;
	
}

}



string toStringLikeInCSharp(double value) {
  import std.format : format;
  return format("%.15G", value);
}
