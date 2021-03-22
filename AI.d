import field;
import analysisEngine;

import std.stdio;
import std.file ;
import std.conv;
import std.getopt;
import std.range;
import std.string;
import std.math;


/* describe the fileformat */
int latPosition 	= 2;
int lonPosition 	= 3;
int timPosition 	= 5;																// position of the element describing time
int stcPosition 	= 4;																// position of the element describing sat count
																						// all positions are given as natural counts,
																						// count does NOT begin at 0 \label countDesc
/* ======================= */

/* describe the input arguments */
string fileName;
int latPosition_updated;
int lonPosition_updated;
int timPosition_updated;
int stcPosition_updated;

/* ======================= */

string category;
double workWidth = -1;
double maxDeviation_ofAngle = -9999;
double minLength = -9999;
double maxDist   = -9999;
double cDist     = -9999;
double geoFenceRadius = -9999;
double [] geoFenceCoordinates = new double [] (0);

void main(string[] args) {

	auto opts = getopt(args, 
						"filename", &fileName,
						"latPos"  , &latPosition_updated, 
						"lonPos"  , &lonPosition_updated,
						"timPos"  , &timPosition_updated,
						"stcPos"  , &stcPosition_updated
								);														
								
	
	if (latPosition_updated != 0 ) latPosition = latPosition_updated;					// see line \ref countDesc. 
																						// because counting does NOT start as zero, 
																						// we can take zero as an indicator 
																						// of missing arguments
																						// we will never say some position = 0, 
																						// because count starts at 1
	if (lonPosition_updated != 0 ) lonPosition = lonPosition_updated;
	if (timPosition_updated != 0 ) timPosition = timPosition_updated;
	if (stcPosition_updated != 0 ) stcPosition = stcPosition_updated;
																						// IN an emergency, we can re_designate 
																						// a different file format
	
	
	field.greenLine [] gLines ;
	field.rawData[] rData ;
	
																						writeln("Extracting values from :" ~fileName);	
	rData = get_rawData_fromFile(fileName);												writeln("Received data points : " ~ to!string(rData.length));
																						// since the file has sequencial waypoints, 
																						// we can guarantee, that rData[n] is a point,
																						// where the vehicle arrived before rData[n+1]
	
	
	
	analysisEngine.geoEngine engine		= new analysisEngine.geoEngine();				// instance began
	
	engine.dataSet				= &rData;												// this data is an array of structs.
																						// elements are lat, lon,satcount, 
																						// and then a string with date and time, 
																						// this one is not set by a set function, 
																						// because of may back and forth casting, 
																						// and back and forth access. 
																						
	engine.set_dataType(2);																// data is array of lats and lons, 
																						// the data type is a private member
	engine.fileName= fileName;
	if(workWidth != -1) engine.workWidth = workWidth;
	if(maxDeviation_ofAngle != -9999) engine.maxDeviation_ofAngle = maxDeviation_ofAngle *PI / 180.00;
	if(minLength != -9999) engine.minLineLength = minLength;
	if(maxDist != -9999) engine.maxLineDist = maxDist;
	if(cDist != -9999) engine.crad = cDist;
	if(geoFenceRadius != -9999) engine.geoFenceRadiusThreshhold = geoFenceRadius;
	if(geoFenceCoordinates.length != 0) engine.geoFenceCoordinates = geoFenceCoordinates;
	
	engine.calculate_drawingParams();
	
	engine.create_map_fromRawData_inDLIB (-1,-1);										// -1, -1 = automatic. 

	//engine.overlay_convexHull_onOutputMap();											// get some idea of where we are going ...
																						writeln("category is : " ~ category);

	final switch (category) {
		
		case "transport": 
			//engine.analyze_transport();
		break;
		
		case "tillage" :																//writeln("skipping");
			engine.analyze_tillage();
		break;
		
		case "grasslandCultivation" :
			engine.analyze_cultivation() ;
		break;
		
	}
	
	//engine.submit_results();
	engine.draw_map();
	
	
	

}





field.rawData[]  get_rawData_fromFile(string fpath) {

	int i = 0;
	
	field.rawData[] allData = new field.rawData[0];

	File f = File(fpath, "r");
	int [] ids = new int[0];
	
	while(!f.eof()) {
		
		field.rawData rdata ;
		
		string l = chomp(f.readln());
		auto elems = split(l , ",");													//writeln(elems); readln();
																						//writeln(category);
		
		if (i == 0) {
			category = elems[0];	
			i++;
		}
		if(elems.length < 2) continue;
		if(elems[0].strip() == "machineWidth") workWidth   = to!double(elems[1].strip());
		if(elems[0].strip() == "slp") maxDeviation_ofAngle  = to!double(elems[1].strip()); 
		if(elems[0].strip() == "lng") minLength  = to!double(elems[1].strip()); 
		if(elems[0].strip() == "dist") maxDist  = to!double(elems[1].strip()); 
		if(elems[0].strip() == "cdist") cDist  = to!double(elems[1].strip()); 
		if(elems[0].strip() == "geoFenceRadius") geoFenceRadius  = to!double(elems[1].strip()); 
		if(elems[0].strip() == "geoFenceCoordinates") geoFenceCoordinates  = [ to!double(elems[1].strip()) ,to!double(elems[2].strip())] ; 
		
		
		if (elems.length < 6) 	continue;
		
		int k = to!int(elems[0].strip());
		
		bool cont = false;
		
		foreach(ii ; ids) {
			if (k == ii) {
				cont = true;
				break;																	
			}
		}
																					// write(k);writeln( " will not be picked up");
		if (cont) continue;
		
		ids ~= k;
		
		rdata.lat = to!double(elems[latPosition -1][0..2]) + to!double(elems[latPosition -1][2..$])/60;	//writeln(rdata.lat);
		rdata.lon = to!double(elems[lonPosition -1][0..3]) + to!double(elems[lonPosition -1][3..$])/60;	//writeln(rdata.lon);
		rdata.satcount = to!short(elems[stcPosition-1]);
		rdata.datestr = elems[timPosition -1];									
		rdata.id = elems[0].strip();
		allData ~= [rdata];
		
		i++;
	}																			
	
	return allData;
	
}
