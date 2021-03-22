module neuron;



import std.stdio;

struct edgeDetectorNeuron {

	double angle 	= 0;
	double width 	= 0;
	long pointCount = 0;
	
	edgeDetectorNeuron_area stripe;
	
	double score 	= 0;
	
	double getScore() {
	
		double s = 0;
		
		return s;
	}
}

struct edgeDetectorNeuron_area {
	double [2] v1	;
	double [2] v2	;
	double [2] v3	;
	double [2] v4	;
	
	int[] detectedPoints;

}
