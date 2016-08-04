library(nlme);

pulseDecayInit = function(mCall, LHS, data){
	# This function estimates good starting parameters for the nls command

	xy=sortedXyData(mCall[["x"]], LHS, data);

	# Greatest Y value
	mag=max(xy[,"y"]);
	# Time of greatest Y value 
	peakt=mean(xy[(which(xy[,"y"]==mag)), "x"]);
	# Standard deviation of time points
	peakw=log(sd(xy[,"x"]));

	value=c(mag, peakt, peakw);
	names(value)=mCall[c("magnitude", "peak_time", "peak_width")];
	return(value);
}

pulseDecayFormula=~ magnitude*exp(-(log(peak_time)-log(x))^2/log(peak_width)^2);

pulseDecayFunction=function(x, magnitude, peak_time, peak_width){
	val= magnitude*exp(-(log(peak_time)-log(x))^2/log(peak_width)^2);
	return(val);
}

#plot(t,pulseDecayFunction(t,3,3,4));

# Define a self start nonlinear function for nls
SSpulseDecay=selfStart(pulseDecayFormula, initial=pulseDecayInit, parameters=c("magnitude", "peak_time", "peak_width"));

if(0){
	x=c(1,5,8,20,33);
	y=c(3,2,3,2,0);
	data_val=cbind(x,y);
	colnames(data_val)=c("testx", "testy");
	df=as.data.frame(data_val);

	getInitial(testy~SSpulseDecay(testx, magnitude, peak_time, peak_width), data=df);
	plot(df);points(t, pulseDecayFunction(t, 5,5, 2.56), type="l");

	nls(testy~SSpulseDecay(testx, magnitude, peak_time, peak_width), data=df)
}
