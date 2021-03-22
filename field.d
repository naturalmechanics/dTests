module field;



import std.stdio;


struct greenLine {

	double startLat 	= 0;
	double startLon 	= 0;
	double endLat		= 0;
	double endLon		= 0;
	double globalSlope	= 0;
	double length		= 0;
	double slopeDelta	= 0;
	double lengthDelta	= 0;
	
	greenLine[][][] parallelLines_withinRange;
	
};


struct rawData {

	double lat		= 0;
	double lon		= 0;
	short satcount		= 0;
	bool 	sensor1		= false, 
			sensor2		= false;
	string datestr		= "";
	string id           = "";

};

struct voronoiRegion {

	double centerLat    = 0;
	double centerLon    = 0;

	int dummyIdx1       = 0;
	int dummyIdx2       = 0;

	double semiMajorAxis= 0;
	double semiMinorAxis= 0;
	double eccentricity = 0;

	double membercount  = 0;
	double relativeFrequency= 0;

	int [] indices_fromRawData = [];
	geoPoint [] actualPoints   = [];
    
    
} ;


struct geoPoint {

	double latitude = 0;
	double longitude= 0;
	

}

struct line {

	long uniqueID       = 0;
	double startLat 	= 0;
	double startLon 	= 0;
	double endLat		= 0;
	double endLon		= 0;
	double slope		= 0;
	double length		= 0;
	ulong startidx 		= 0;
	ulong endidx		= 0;
    
}

class field {

	this(){}
	~this(){}

    line[] borderLines;
    line[] cultivatedLines;
    line[] ignoredLines;
	
}

class cluster {

}

