#include <math.h>

#include "./GravityTensor/BaseDefine.h"
#include "./GravityTensor/GravityVW.h"

int main()
{

	double R[3];
	double lon,lat,alt;
	double U,Ux,Uy,Uz,Uxx,Uyy,Uzz,Uxy,Uxz,Uyz;
	int Flag;

	CGravityTensor  GTensor;

	GTensor.SetGrData(&GrData_WGS84_EGM96,&FullNormalFactor);
	GTensor.SetDegreeOrder(2,4,4);

	lon = 0;
	lat = 90;
	alt = 7378140;

	while(lat>=-90)
	{
		R[0] = alt*cos(lat*RADDEG)*cos(lon*RADDEG);
		R[1] = alt*cos(lat*RADDEG)*sin(lon*RADDEG);
		R[2] = alt*sin(lat*RADDEG);

		GTensor.SetPosition(R);

		Flag = GTensor.GetTensorValue(0,0,0,&U);
		Flag = GTensor.GetTensorValue(1,0,0,&Ux);
		Flag = GTensor.GetTensorValue(0,1,0,&Uy);
		Flag = GTensor.GetTensorValue(0,0,1,&Uz);
		Flag = GTensor.GetTensorValue(2,0,0,&Uxx);
		Flag = GTensor.GetTensorValue(0,2,0,&Uyy);
		Flag = GTensor.GetTensorValue(0,0,2,&Uzz);
		Flag = GTensor.GetTensorValue(1,1,0,&Uxy);
		Flag = GTensor.GetTensorValue(1,0,1,&Uxz);
		Flag = GTensor.GetTensorValue(0,1,1,&Uyz);

		double a = Uxx+Uyy+Uzz;
		lat = lat-1;
	}

	return 1;

}
