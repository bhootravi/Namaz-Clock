#include<stdio.h>
#include<limits.h>
#include<math.h>
#include <avr/interrupt.h>

#define convention Karachi
#define asrMethod shafai
#define maghribMethod mInterval

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

#define sunriseInterval (10/60.0)
#define ishraqInterval (15/60.0)
#define zawalInterval (5/60.0)
#define sunsetInterval (3/60.0)
#define sahriInterval (15/60.0)
#define iftarInterval (15/60.0)

//shafai includes Shafi'i, Maliki, Ja'fari, and Hanbali
#define shafai 1
#define hanafi 2

#define mInterval 1
#define mDirect 2

enum events{
			fajr, sunrise, ishraq, chast, zawal,
			zuhar, asr, sunset,
			maghrib, isha, tahajjud, sahri, iftar
} event;

float angles[7][2] = {{18,17},{15,15},{19.5,17.5},{18.5,-1},{18,18},{17.7,14},{16,14}};

double start[13];
int hm[13][2];
int hijri_date[3];

int time[3], greg_date[3];
unsigned char days[13] = {-1,31,28,31,30,31,30,31,31,30,31,30,31};

// double end[13];

#define set_flag(bit) (flags |= (1 << (bit)))
#define check_flag(bit) (flags & (1 << (bit)))
#define clear_flag(bit) (flags &= !(1 << (bit)))

//bit 0 - MIN_changed
#define MIN_changed 0
//bit 1 - HOUR_changed
#define HOUR_changed 1
//bit 2 - DATE_changed
#define DATE_changed 2
//bit 3 - LEAP_YEAR
#define LEAP_YEAR 3
unsigned char flags = 0;

double EqT, declination;

// double 

ISR(TIMER1_COMPA_vect)
{
	time[2]++;
	if(time[2] == 60)
	{
		time[2] = 0;
		set_flag(MIN_changed);
	}
}

void update_time(void)
{
	time[1]++;
	if(time[1] == 60)
	{
		time[1] = 0;
		time[0]++;
		set_flag(HOUR_changed);
	}
	if(time[0] == 24)
	{
		time[0] = 0;
		greg_date[0]++;
		set_flag(DATE_changed);
	}
	unsigned char t = (check_flag(LEAP_YEAR) && greg_date[1] == 2) ? 1:0;
	if(greg_date[0] > (days[greg_date[1]]+t))
	{
		//TODO february leap year
		greg_date[0] = 1;
		greg_date[1]++;
	}
	if(greg_date[1] == 13)
	{
		greg_date[1] = 1;
		greg_date[2]++;
		if((greg_date[2] & 3) == 0 && ((greg_date[2] % 25) != 0 || (greg_date[2] & 15) == 0))
			set_flag(LEAP_YEAR);
	}
}

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
	return ((long)(2 - a + (a / 4)) + d + (long)(365.25 * (y + 4716)) 
	+ (long)(30.6001 * (m + 1)) - 1524);
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
	
	event = ishraq;
	start[event] = start[sunrise] + ishraqInterval;
	
	event = chast;
	start[event] = (start[sunrise] + start[zawal]) / 2;
	
	event = zawal;
	start[event] = noonT;
	dtohm(start[event], hm[event]);
	
	event = zuhar;
	start[event] = noonT + zawalInterval;
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
	if(maghribMethod == mInterval)
	{
		start[event] = start[sunset] + sunsetInterval;
	}
	else if(maghribMethod == mDirect)
	{
		start[event] = noonT + sunAngle(4);
	}
	dtohm(start[event], hm[event]);
	
	event = isha;
	if(angles[convention][1] > 0)
	{
		start[event] = noonT + sunAngle(angles[convention][1]);
		dtohm(start[event], hm[event]);
	}
	else
	{
		//TODO 120 mins during Ramadan
		start[event] = start[maghrib] + 1.5;
	}
	
	event = tahajjud;
	start[event] = 0.0;
	
}

void hijri(int date, int month, int year, int* out)
{
	unsigned long jd = julian(date, month, year);
	
	long l = jd - 1948440 + 10632;
	long n = (long)((l - 1) / 10631.0);
	l = l - 10631 * n + 354;
	long j = ((long)((10985 - l) / 5316.0)) * ((long)((50 * l) / 17719.0));
	j += ((long)(l / 5670.0)) * ((long)((43 * l) / 15238.0));
	long l1 = ((long)((30 - j) / 15.0)) * ((long)((17719 * j) / 50.0));
	l1 +=  ((long)(j / 16.0)) * ((long)((15238 * j) / 43.0));
	l = l - (l1 - 29);
	out[1] = (long)((24 * l) / 709.0);
	out[0] = l - (long)((709 * out[1]) / 24.0);
	out[2] = 30 * n + j - 30;
}

void hijri_2(int date, int month, int year, int* out)
{
	int dayear[13] = {-1, 31, 28, 31, 30, 31, 30 , 31, 31, 30, 31, 30, 31};
	int mday=0;
	int i;
	for(i = 1; i < month; i++)
		mday = mday + dayear[i];
	if((year % 4) == 0)
	{
		if(year % 100 == 0)
		{
			if(year % 400 == 0)
				mday++;
		}
		else
			mday++;
	}
	double greg = year + ((mday + date) / 365.0);
	double islm = ((greg - 621.5774)*1000000)/ 970224.0;
	
	int islmyear = (int)(islm);
	int islmday = (int)((islm - islmyear) * 360);
	
	int ayear = islmyear;
	int aday = islmday;	
	int amonth = 0;
	
	while (aday > 0)
	{
		amonth = amonth + 1;
		aday = aday - 30;
	}
	aday = aday + 33;
	if (aday > 30)
	{
		aday = aday - 30;
		amonth = amonth + 1;
	}
	if (amonth == 13)
	{
		amonth = 1;
		ayear = ayear + 1;
	}
	out[2] = ayear;
	out[1] = amonth;
	out[0] = aday;
}

void RTC_init(void)
{
	//Do this first
	//Warning: When switching between asynchronous and synchronous clocking of
	//Timer/Counter2, the Timer Registers TCNT2, OCR2, and TCCR2 might be corrupted.
	//External Crystal
	ASSR |= 0x08;
	// Prescalar = f/1024
	// Enable CTC mode
	TCCR2 |= 0x0F;
	// 2^15 ticks = 1 s for 32.768 kHz crystal
	// Prescalar f/1024 = 2^10 ticks
	OCR2 = 32;
	//Wait for the settings to be updated
	while((ASSR & 0x3) != 3);
	//Enable Output Compare Match Interrupt 
	TIMSK |= 0x80;
}

int main(int argc, char* argv[])
{
	greg_date[0] = 26;
	greg_date[1] = 8;
	greg_date[2] = 2015;
	int date = greg_date[0], month = greg_date[1], year = greg_date[3];
	int h = time[0], m = time[1], s = time[2];
	
	RTC_init();
	sei();
	if(check_flag(MIN_changed))
	{
		update_time();
		update_event();
		clear_flag(MIN_changed);
	}
	if(check_flag(DATE_changed))
	{
		calcTimes(date, month, year);
		//to match calender from http://www.makkahcalendar.org/
		if((greg_date[1]%2) == 0)
			hijri(date, month, year, hijri_date);
		else
			hijri_2(date, month, year, hijri_date);
		clear_flag(DATE_changed);
	}

// 	int zh = (int)zawal;
// 	int zm = (int) ((zawal - zh) * 60);
	return 0;
}

