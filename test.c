// Sunrise and Noon OK

#include<stdio.h>
#include<limits.h>
#include<math.h>

#define convention Karachi
#define asrMethod shafai
double latitude = 22.2833;
double longitude = 70.7833;
//in m
int elevation = 300;
float timezone = 5.5;

#define radians(d) (d)*M_PI/180.0
#define degrees(r) (r)*M_1_PI*180.0

#define MWL 0
#define ISNA 1
#define Egypt 2
#define Makkah 3
#define Karachi 4
#define Tehran 5
#define Qum 6

#define ishraqInterval 15

//shafai includes Shafi'i, Maliki, Ja'fari, and Hanbali
#define shafai 1
#define hanafi 2

enum events{
			fajr, sunrise, ishraq, chast, zawal,
			zuhar, asr, sunset,
			maghrib, isha, tahajjud, sahri, iftar
} event;

float angles[7][2] = {{18,17},{15,15},{19.5,17.5},{18.5,-1},{18,18},{17.7,14},{16,14}};

int hm[13][2];
double start[13];
double end[13];

double EqT, declination;



// double 

double normalise(double n, double m)
{
	return n - (m * (int)(n / m));
}

double acot(const double x)
{
	return atan2(1.0, (x));
}

// source http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
unsigned long julian(int d, int m, int y)
{
	if(m <= 2)
	{
		y -= 1;
		m += 12;
	}
	int a = y/100;
	//for integral ans(noon) -1524
 	return ((int)(2 - a + (a / 4)) + d + (int)(365.25 * (y + 4716)) 
 			+ (int)(30.6001 * (m + 1)) - 1524);
}

void solarPosition(int date, int month, int year)
{
	int d = julian(date, month, year) - 2451545;
	//Mean anomaly of the Sun
	double g = 357.529 + 0.98560028 * d;
	g = normalise(g, 360);
	//Mean longitude of the Sun:
	double q = 280.459 + 0.98564736 * d;
	q = normalise(q, 360);
	//Geocentric apparent ecliptic longitude of Sun(adjusted for aberration)
	double L = q + 1.915 * sin(radians(g)) + 0.020 * sin(2 * g);
	L = normalise(L, 360);
	//The distance of the Sun from the Earth, R, in AU (approx.)
// 	double R = 1.00014 - 0.01671 * cos(radians(g)) - 0.00014 * cos(radians(2 * g));
	//mean obliquity of the ecliptic, in degrees
	double e = 23.439 - 0.00000036 * d;
	//Sun's right ascension
	double RA=degrees(atan2(cos(radians(e))*sin(radians(L)),cos(radians(L))));
	//convert to hours
	RA /= 15.0;
	RA = normalise(RA, 24);
	//Sun's declination
	declination = degrees(asin(sin(radians(e)) * sin(radians(L))));
	//The Equation of Time, EqT, apparent solar time minus mean solar time
	EqT = q/15.0 - RA;
	//The angular semidiameter of the Sun, SD, in degrees
// 	double SD = 0.2666 / R;
}

//a - angle of sun with horizon, return time difference from noon
double sunAngle(float a)
{
	//source http://praytimes.org/calculation/
	if(elevation > 0)
		a += 0.0347 * sqrt(elevation);
	double n = (sin(radians(latitude)) * sin(radians(declination)));
	n = -sin(radians(a)) - n;
	double d = cos(radians(latitude)) * cos(radians(declination));
	double t = degrees(acos(n / d));
	t /= 15;
	return t;
}

double noon(void)
{
	//http://praytimes.org/calculation/
	return (12 + timezone - longitude/15 - EqT);
}

// return hours and mins from decimal repr.
void dtohm(double t, int* hm)
{
	hm[0] = (int) t;
	hm[1] = (int) ((t - hm[0]) * 60);
}

void calcTimes(int date, int month, int year)
{
	solarPosition(date, month, year);
	double noonT = noon();
	event = fajr;
	start[event] = noonT - sunAngle(angles[convention][0]);
	dtohm(start[event], hm[event]);
	event = sunrise;
	start[event] = noonT - sunAngle(0.8333);
	dtohm(start[event], hm[event]);
	event = zawal;
	start[event] = noonT;
	dtohm(start[event], hm[event]);
	event = asr;
	double temp = acot((asrMethod + tan(latitude - declination)));
	start[event] = noonT + sunAngle(-temp);
	dtohm(start[event], hm[event]);
	event = sunset;
	start[event] = noonT + sunAngle(0.8333);
	dtohm(start[event], hm[event]);
	//other method 
	//sunset time + 1-2 mins.
	event = maghrib;
	start[event] = noonT + sunAngle(4);
	dtohm(start[event], hm[event]);
	event = isha;
	if(angles[convention][1] > 0)
	{
		start[event] = noonT + sunAngle(angles[convention][1]);
		dtohm(start[event], hm[event]);
	}
	//TODO else for Makkah
}

int main(int argc, char* argv[])
{
	int date = 14, month = 7, year = 2015;
	calcTimes(date, month, year);

// 	int zh = (int)zawal;
// 	int zm = (int) ((zawal - zh) * 60);
	return 0;
}

