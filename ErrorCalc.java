import java.lang.Math;
import java.util.Arrays;
public class ErrorCalc{
	public static void main(String[] args)
	{
		int length = 1000;
		double B = 4131; //4131K but this is between temps of 25 Cel and 100 Cel not the min and max temp
		double RTemp = 47; //R25 resistance at 25 degree cel.. THIS IS 45 KOHMS
		double minTemp = 273; // should be in Kelvin
		double maxTemp = 300; //
		double VCC = 5.5; //Assuming a 5.5 supply voltage?
		double ROne = 47; //resistor divider..resistance is 47kOhms
		double VoutMax = Vout(maxTemp, B, RTemp, VCC, ROne);
		double VoutMin = Vout(minTemp, B, RTemp, VCC, ROne);
		double Range = Math.abs(VoutMax-VoutMin);
		System.out.println("Range is " + Range);
		System.out.println("VoutMax is " + VoutMax);
		System.out.println("VoutMin is " + VoutMin);
		double[] RTArray = new double[length];
		double[] TArray = new double[length];
		double[] EThermArray = new double[length];
		double[] EresArray = new double[length];
		double[] EADCArray = new double[length];
		double[] EVCCArray = new double[length];
		double[] VoutIn = new double[length];
		double Delta = Range/length;
		double RTPrime = 0;
		double TError = 0;
		double[] worstE = new double[16]; //worst error
		double MaxError = 0;
		double indexWC;
		double Increase= VoutMin;
		double RTol = 0.5; //0.5%
		double BTol = .5; //0.5% WHAT TO SET?
		double R = 47; //for ERES the resistor value
		double Etemp = 10*Math.pow(10,-6); //the temperature coefficent [ppm/Celcius]
		double ETol = .05; //the resistance for Resistor Tolerance
		double INL = 10; //INL, FSR, EGain, EOff, EOffD, EGainD, are from the data sheet
		double FSR = 4.096;
		double EGain = 0; //not sure.. says we can select x1, x2...
		double EOff = 30;
		double EOffD = 300;
		double n = 16; //number of bits in ADC... how to find?
		double EQ = 0; //FSR/(Math.pow(2,n)); //use equation to get this... n is the number of bits in ADC
		double EGainD = 15;
		double percentageError = .04;
		double EOutD = 10*Math.pow(10,-6); //cannot find
		double ELoad = (5.5*Math.pow(10,-6))/50;
		double[] DiffTemp = new double[1000];
		double MaxDiffTemp;
		System.out.println("Delta " + Delta);
		System.out.println("Increase " + Increase);
		System.out.println("Increase/VCC " + Increase/VCC);

//ask what type of resistor we are using.. ERA 1A? how to know which temperature coefficent
		for(int i=0; i<length; i++)
		{
			RTArray[i]= RT(Increase, ROne, VCC);
			TArray[i] = TempCalc(RTArray[i], RTemp, B);
			//TArray[i] = 25; //For Debugging set TArray elements to 25
			EThermArray[i] = ETherm(RTemp, RTol, B, BTol, TArray[i]);
			EresArray[i] = Eres(R, Etemp, TArray[i], ETol);
			EADCArray[i] = EADC(INL, FSR, EGain, EOff, EOffD, EQ, EGainD, TArray[i]);
			EVCCArray[i] = EVCC(percentageError, VCC, TArray[i], RTArray[i], R, EOutD, ELoad);
			VoutIn[i] = Increase;
			Increase = Delta + Increase; //because we are incrementing between the ranges of VOUT
		}

		for(int i=0; i<16; i++)
		{
			for(int j=0; j<1000; j++)
			{

		//Equation 12, R'T we found then plus Etherm
		byte BitOne = (byte)0b00000001; //1st two loops first post then neg to repeat for each of the bits
		byte BitTwo = (byte)0b00000010;
		byte BitThree = (byte)0b00000100;
		byte BitFour = (byte)0b00001000;

		RTPrime = RT(VoutIn[j] + (EADCArray[j]/(1*Math.pow(10,6)))*(1-(2*(BitOne&i))),ROne + EresArray[j]*(1-(2*(BitTwo&i) >> 1)), VCC + EVCCArray[j]*(1-(2*(BitThree&i) >> 2)));
		TError = TempCalc(RTPrime+EThermArray[j]*(1-(2*(BitFour&i) >> 3)), RTemp, B);

		DiffTemp[j] = Math.abs(TArray[j] - TError); //Array of The Change in T


		  }

			System.out.println("EVCCArray " + EVCCArray[length-1]);

			System.out.println("EADCArray " + EADCArray[length -1]);

			System.out.println("RT PRIME " + RTPrime);

			System.out.println("TError " + TError);

			System.out.println("EThermArray " + EThermArray[length-1]);

			System.out.println("TArray " + TArray[length-1]);

			System.out.println("RTArray " + RTArray[length-1]);

			System.out.println("ERESArray " + EresArray[length-1]);

			MaxDiffTemp = Max(DiffTemp);
				if(MaxDiffTemp >= MaxError) //max of DiffTEmp[i] replaces diffTemp
				{
						worstE = DiffTemp; //would it be worstE[i] = DiffTemp[i] or worstE = DiffTemp[i]
						MaxError = MaxDiffTemp;
						indexWC = i;
				}

		}

		System.out.println("Worst Error " + Arrays.toString(worstE)); //in degrees
	}

	public static double Max(double array[])
	{
		double max = 0;
		for(int i =0; i<array.length; i++)
		{
			if(array[i] > max)
			{
				max = array[i];
			}
		}
		return max;
	}

	public static double Vout(double temp, double B, double RTemp, double VCC, double ROne) //equation 1// Min + Max
	{
		double RT = Bmethod(temp, B, RTemp);
		double Vout = VCC*(ROne/(ROne+RT));
		return Vout;
	}

	public static double RT(double MeasVout, double ROne, double VCC)
	{

		/*
		double a = MeasVout/VCC;
		double RT = (ROne*(1-a))/a;
		*/

		//double RT = ROne*(VCC/MeasVout) - ROne;

		double RT = ((ROne*VCC)/MeasVout) - ROne;
		return RT; //is this wrong
	}

	public static double Bmethod(double Temperature, double B, double RTemp) //Min+Max
	{
		double RT;
		double K = 298.15; //Taking our temp outside and converting to Kelvins
		RT = RTemp * Math.exp(B*((1/Temperature) - (1/K))); //RTemp = R25
		//Ask Why we have an RT equation for finding Resistance
		//Ask why we need T equation for using the resistance we found to find temperature
		return RT;

	}

	public static double TempCalc(double RT, double RTemp, double B)
	{
		double K = 298.15;
		double T = B/(Math.log(RT/RTemp)+(B/K)); //naturla log
		return T;
	}

	public static double Eres(double R, double Etemp, double Temp, double ETol)
	{
		//Etol and E temp are constants
		//Consider it absolute value
		double Eres = ((ETol/100)*R)+(Etemp*R*Temp); //do not need 1e6
		return Eres;
	}

	public static double EVCC(double percentageError, double VCC, double RTemp, double RT, double R, double EOutD, double ELoad)
	{
		double EVCC = (percentageError*VCC) + (EOutD*VCC) + ELoad*(VCC/(RTemp+R));
		return EVCC;
	}

	public static double EADC(double INL, double FSR, double EGain, double EOff, double EOffD, double EQ, double EGainD, double T)
	{
		double EADC = Math.sqrt(Math.pow((INL*FSR),2) + Math.pow(((EGain/100)*FSR+EOff),2) + Math.pow(((EOffD/1000)*T),2)+ Math.pow((EGainD*T),2)+Math.pow(EQ,2));
		return EADC;
	}

	public static double ETherm(double RTemp, double RTol, double B, double BTol, double temp)
	{
		//RTolTherm same name but different number
		double ETherm = Bmethod(temp,(B*BTol)/100, (RTol*RTemp)/100);
		return ETherm;
	}
}
