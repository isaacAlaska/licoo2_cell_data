/**
  C language lipo battery simulator
  
  Written by Dr. Orion Lawlor <lawlor@alaska.edu> 2018-03-07 (Public Domain)
  Approach and calibration parameters from Isaac Thompson's MS Thesis 2018.
*/
#include <stdio.h>

/*
Battery model for rechargable lithium-ion cell.
*/
struct battery_model {
	/* Fully charged capacity, in amp-seconds */
	float capacityAs;
	
	/* State of charge, from 0.0 (fully discharged) to 1.0 (fully charged) */
	float SOC;
	
	/* Charge (coloumbs) borrowed from short term capacitor C1 */
	float C1Q;
	
	/* Temperature (deg C) of interior of cells */
	float cellT;
};



/* Stores battery model parameters for the current battery configuration */
struct battery_model_parameters {
	float Em; /* Open circuit voltage (volts) */
	float R0; /* Series output resistance (ohms) */
	float R1; /* Short term deep draw resistance (ohms) */
	float C1; /* Short term capacitance (farads) */
};


/* Stores a table of one parameter for the battery model 
   where the table is indexed by battery SOC and temperature. */
struct battery_model_table {
#define battery_model_table_SOCs 11  /* number of entries in table by state of charge: 0.0 0.1 ... 1.0 */
#define battery_model_table_temps 4 /* number of entries in table by temperature */
	float values[battery_model_table_temps][battery_model_table_SOCs];
};

/* Temperature values */
const static float battery_model_temperatures[battery_model_table_temps]={-20.0, -10.0, -5.0, +2.0};

/* Open circuit voltage, Em (volts) */
const static struct battery_model_table battery_model_Em={{
	{3.5,3.65,3.7,3.75,3.78,3.8,3.85,3.9,3.95,4.1,4.2}, // -20 deg C
	{3.5,3.65,3.7,3.746368,3.794009,3.824597,3.870755,3.921037,3.984153,4.1,4.2}, // -10 deg C
	{3.5,3.717802,3.751656,3.779548,3.805342,3.837747,3.886275,3.92452,4.019383,4.131402,4.2}, // -5 deg C
	{3.5,3.723299,3.754516,3.788628,3.812054,3.840599,3.888213,3.933897,4.024288,4.130746,4.182739}, // 2 deg C
}};
/* Series output resistance R0 (ohms) */
const static struct battery_model_table battery_model_R0={{
	{0.26,0.26,0.26,0.13,0.13,0.13,0.13,0.13,0.25,0.2,0.67}, // -20 deg C
	{0.3,0.050589,0.144401,0.085073,0.091675,0.085872,0.08382,0.084737,0.075961,0.15,0.25}, // -10 deg C
	{0.2,0.029142,0.029737,0.031219,0.031587,0.030885,0.031477,0.030845,0.030875,0.025,0.016}, // -5 deg C
	{0.032564,0.022225,0.019854,0.024638,0.022878,0.021342,0.022003,0.02195,0.021421,0.023454,0.014168}, // 2 deg C
}};
/* Short term deep draw resistance R1 (ohms) */
const static struct battery_model_table battery_model_R1={{
	{2,0.75,0.21,0.190953,0.147748,0.127334,0.143009,0.180778,0.1,0.261743,0.85}, // -20 deg C
	{0.003815,0.007988,0.020238,0.015108,0.01404,0.014878,0.014838,0.014781,0.015083,0.15,0.3}, // -10 deg C
	{0.011421,0.003253,0.012514,0.00939,0.010378,0.009284,0.008821,0.008391,0.010644,0.008414,0.007233}, // -5 deg C
	{0.025991,0.003294,0.013872,0.013772,0.013957,0.011306,0.01088,0.01135,0.015937,0.012274,0.007585}, // 2 deg C
}};
/* Short term capacitance C1 (farads) */
const static struct battery_model_table battery_model_C1={{
	{400,500,600,846,846,846,846,846,600,846,596}, // -20 deg C
	{14.34898,28719.38,1818.858,5832.355,8962.667,8772.705,8750.688,8565.881,7004.807,11188.4,7370.326}, // -10 deg C
	{0.881527,33414.97,2179.029,11289.18,7234.158,6226.428,5750.18,9030.291,3869.932,11851,7122.03}, // -5 deg C
	{0.262732,50759.86,3022.06,15720.72,8308.124,7180.572,6619.685,13150.94,4201.662,15103.12,6852.036}, // 2 deg C
}};

/* Bilinear interpolation of one parameter from this table of battery parameters */
float battery_model_interpolate(const struct battery_model_table *table,
	float T_number,int T_index,float SOC_number,int SOC_index)
{
	int SOC_next=SOC_index+1;
	if (SOC_next>=battery_model_table_SOCs) SOC_next=battery_model_table_SOCs-1;
	int T_next=T_index+1;
	if (T_next>=battery_model_table_temps) T_next=battery_model_table_temps;
	float II=table->values[T_index][SOC_index];
	float IN=table->values[T_index][SOC_next];
	float TI=table->values[T_next ][SOC_index];
	float TN=table->values[T_next ][SOC_next];
	float I=II + (IN-II)*(SOC_number-SOC_index);
	float T=TI + (TN-TI)*(SOC_number-SOC_index);
	float ret=I + (T-I)*(T_number-T_index);
	// printf(" Interpolating T=%f S=%f: %f\n",T_number,SOC_number,ret);
	return ret;
}

/* Look up the currently applicable entry in this battery model parameter table. */
void battery_model_get_parameters(const struct battery_model *battery,struct battery_model_parameters *param)
{
	// State of charge table is distributed uniformly
	float SOC_number=(battery->SOC)*(battery_model_table_SOCs-1);
	int SOC_index=(int)SOC_number;
	if (SOC_index<0) { SOC_number=SOC_index=0; }
	if (SOC_index>=battery_model_table_SOCs) { SOC_number=SOC_index=battery_model_table_SOCs-1; }
	
	// Look up temperature in table of temperatures:
	float lookupT=battery->cellT;
	int T_index=0;
	while (T_index+1<battery_model_table_temps 
	    && battery_model_temperatures[T_index+1]<=lookupT)
			T_index++;
	float T_number=0;
	if (T_index+1<battery_model_table_temps) {
		// linearly interpolate between nearest temperatures
		float last=battery_model_temperatures[T_index];
		float next=battery_model_temperatures[T_index+1];
		T_number=T_index + (lookupT-last)/(next-last);
	}
	
	param->Em=battery_model_interpolate(&battery_model_Em,T_number,T_index,SOC_number,SOC_index);
	param->R0=battery_model_interpolate(&battery_model_R0,T_number,T_index,SOC_number,SOC_index);
	param->R1=battery_model_interpolate(&battery_model_R1,T_number,T_index,SOC_number,SOC_index);
	param->C1=battery_model_interpolate(&battery_model_C1,T_number,T_index,SOC_number,SOC_index);
}

/*
 Battery model circuit:
	Idealized voltage source Em
	Parallel short-term resistance R1 and capacitance C1
	Series output resistor R0
*/

/* Create a new battery model with the given:
	Capacity, in amp hours
	state of charge (0.0 to 1.0 fully charged)
	and temperature (deg C) 
*/
void battery_model_init(struct battery_model *battery,float capacityAh,float SOC,float tempC)
{
	battery->capacityAs=capacityAh*3600.0;
	battery->SOC=SOC;
	battery->cellT=tempC;
	battery->C1Q=0.0; // assume C1 is at equilibrium
}

/* Estimate the voltage (volts) at the battery output terminals
   that the battery will supply at this draw current (amps). */
float battery_model_voltage(const struct battery_model *battery,float amps)
{
	struct battery_model_parameters param;
	battery_model_get_parameters(battery,&param);
	
	/* voltage drop across R0 */
	float R0V=param.R0*amps; 

	/* voltage drop across resistor R1 = voltage drop across capacitor C1 */
	float R1V=battery->C1Q/param.C1;
	
	return param.Em - R1V - R0V;
}


/* Update the battery electrical state based on this current draw for this time. 
      amps is the measured current draw (amperes)
      dt is the simulation timestep (seconds)
   Returns the heat energy, in Joules, added to the battery.
*/
float battery_model_electrical(struct battery_model *battery,float amps, float dt)
{
	struct battery_model_parameters param;
	battery_model_get_parameters(battery,&param);
	
	float R0I=amps;
	float R0V=param.R0*R0I;
	
	float C1V=battery->C1Q/param.C1; // voltage across C1
	float R1V=C1V; // voltage across R1
	float R1I=R1V/param.R1; // current through R1
	float C1I=amps-R1I; // current flowing out of C1
	
	// printf("		C1V: %.2f V\n",C1V);
	
	battery->C1Q += C1I * dt; // coloumbs of charge leaving C1
	float SOC_amps = amps; // measures SOC after C1
	// SOC_amps = R1I; // measures SOC before C1 (doesn't match reality: eliminates voltage rebound)
	battery->SOC -= SOC_amps * dt / battery->capacityAs;  // coloumbs leaving battery

	// Compute heat emitted by the battery's electrial operation over this period
	float power = R0V*R0I + R1V*R1I;
	float energy = power * dt;
	
	return energy;
}

/* Update the battery heating model:
	heating_joules is the electrical heat energy input, from battery_model_electrical (J)
	specific_heat is the battery's specific heat capacity (joules/(degree C * gram))
	mass is the battery's mass (grams)
	
	ambientT is the ambient temperature (degrees C)
	Rvalue is the battery compartment insulation metric R-value (m^2*degrees C/watt)
	area is the area of the battery compartment exposed to ambient (m^2)
	
	dt is the simulation timestep (seconds)
*/
void battery_model_thermal(struct battery_model *battery,
	float heating_joules, float specific_heat, float mass, 
	float ambientT, float Rvalue, float area, 
	float dt)
{
	float cool_joules=(battery->cellT-ambientT) * area / Rvalue * dt;
	float netT=(heating_joules-cool_joules)/(specific_heat*mass);
	battery->cellT += netT;
}


/* Demo of how to call the simulator
  (modeled after Figure 4-8, -20C where self-heating causes voltage to rise under load)
*/
int main() {
	float ambientT=-20.0;
	struct battery_model battery;
	battery_model_init(&battery,1.8, 1.0, ambientT);
	int S=1; // cells stacked in series
	int dt=12.0; // seconds per timestep
	for (float time=0.0;time<30.0*60.0;time+=dt) {
		float amps=1.8;  // discharge current
		float minutes_between=17.0; // minutes between charge cycles
		float minutes_charge=5.0; // minutes to keep charging at each cycle
		float time_cycle=fmod(time,minutes_between*60.0); // 20 minute charge cycle
		if (time_cycle<10.0 || time_cycle>minutes_charge*60.0+10.0) amps=0.0; // outside charge time
		
		float volts=S*battery_model_voltage(&battery,amps);
		float heat=S*battery_model_electrical(&battery,amps,dt);
		battery_model_thermal(&battery,
			heat, 0.9 /* aluminum J/g */, 150.0 /* grams */, 
			ambientT /* degrees C ambient */, 0.1 /* R value, air film */, 0.1*0.1 /* surface area */,
			dt
		);
		
	//	if (amps>0.0 || fmod(time,60.0)<=0.0)
			printf("%.2f minutes: %.2f V @ %.2f A ( %.2f deg C, %.2f SOC, %.0f C1Q)\n",
				time/60.0, volts, amps, battery.cellT, battery.SOC, battery.C1Q);
	}
	return 0;
}
