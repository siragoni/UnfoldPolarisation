#include <math.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include "TView.h"
#include "TCanvas.h"
#include "TPolyMarker.h"
using namespace std;

TCanvas *c;
TPolyMarker **p;

class V
{
	public:
		V(double a=0, double b=0):x(a),y(b){}
		V(const V& v){ x = v.x; y = v.y; }
		//void setx(double a) { x = a; };
		//void sety(double b){ y = b; };
		double getx(){ return x; };
		double gety(){ return y; };
		double modulo(){ return sqrt(x*x+y*y); };
		V normalizzazione(){ V u(x/modulo(), y/modulo()); return u; }
		V operator+(const V &v){ V u(x+v.x, y+v.y); return u; }
		V& operator+=(const V &v){ x += v.x; y += v.y; return *this; }
		V operator-(const V &v){ V u(x-v.x, y-v.y); return u; }
		//V& operator-=(const V &v){ x -= v.x; y -= v.y; return *this; };
		double operator*(const V &v){ return x*v.x + y*v.y; };
		V operator*(const double &c){ V u(x*c, y*c); return u; }
		V operator/(const double &c){ V u(x/c, y/c); return u; }
	protected:
		double x;
		double y;
};

	const double G=6.672e-11;

class crpCeleste
{

	public:
		crpCeleste(double mass=0, V x0=0, V v0=0);
		void setMass(double mass) {M=mass;};
		void setx(V x0) {x=x0;};
		void setv(V v0) {v=v0;};
		V getx() {return x;};
		V getv() {return v;};
		double m() {return M;};
		V Ftotale(vector<crpCeleste*> lista);
		virtual void evolve(double deltaT, V F){V a=F/M; v+=a*deltaT; x+=v*deltaT;}
		void evolve(double deltaT, vector<crpCeleste*> list){evolve(deltaT,Ftotale(list));}

	protected:
		V x;
		V v;
		double M;
};

	crpCeleste::crpCeleste(double mass, V x0, V v0): M(mass), x(x0), v(v0) {}

	V crpCeleste::Ftotale(vector<crpCeleste*> list)
		{
			 V R;
			 for (unsigned int i=0; i<list.size(); i++)
				{
				 V r = list[i]->getx()-x;
				 double d = r.modulo();
				 V versore = r.normalizzazione();
				 if (d!=0) R += versore*G*M*list[i]->m()/(d*d);
				 }
			 return R;
		 }

class SSolare
{
	public:
		SSolare():deltaT(0),TMax(0){corpi.clear();}
		int how_many() {return corpi.size();}
		void add(crpCeleste *corpo){corpi.push_back(corpo);}
		void setMaxTime(double T) {TMax=T;}
		void setDeltaT(double T) {deltaT=T;}
		void evolvi(double deltaT){for(unsigned int i=0;i<corpi.size();i++)corpi[i]->evolve(deltaT, corpi);}
		virtual void evolvi(void (*f)(vector<crpCeleste*>));
	private:
		double TMax;
		double deltaT;
		vector<crpCeleste *> corpi;
};


void SSolare::evolvi(void (*f)(vector<crpCeleste *>))
{
	double T=0;
	while(T<TMax)
		{
		 evolvi(deltaT);
		 (*f)(corpi);
		 T+=deltaT;
		 }
}

void output(vector<crpCeleste *> corpi)
{
	vector<crpCeleste *>::const_iterator corpo = corpi.begin();
	int count = 0;
	while (corpo != corpi.end())
		 {
			double x = (*corpo)->getx().getx()*400/1e15;
			double y = (*corpo)->getx().gety()*400/1e15;
			p[count]->SetPoint(0, x, y);
			float size = 0.2*(log10((*corpo)->m())-20);
			if (size<0) size = .5;
			p[count++]->SetMarkerSize(size);
			corpo++;
			}
	c->Modified();
	c->Update();
}


class Apollo : public crpCeleste
{
	public:
		Apollo(double mass=0, V x0=0, V v0=0) : crpCeleste(mass, x0, v0) {/*o=origin; Tstart=start;*/ init(); }
		void setOrigin(crpCeleste *origin) { o = origin; }
		void setStartTime(double start) { Tstart = start; }
		//void setrazzi(double T, V v) { vrazzi=v;}
		void accendiRazzi(double T, V v) { Taccensioni.push_back(T+Tstart); vRazzi.push_back(v); }
		virtual void evolve(double dT, V F);
	private:
		/*void init() { Tstart = -1; T = 0; starting = 1; control=1; };
		V vrazzi;
		crpCeleste *o;
		double Tstart;
		double T;
		char starting;
		char control;*/
		void init() { Tstart = -1; T = 0; starting = 1; Nstep =0; };
		vector<double> Taccensioni;
		vector<V> vRazzi;
		crpCeleste *o;
		double Tstart;
		double T;
		char starting;
		unsigned int Nstep;
};

void Apollo::evolve(double dT, V F)
{
	T += dT;
	if (Tstart < 0) return;
	if (T < Tstart) x = o->getx();
	else
		{
			if (starting)
				{
			      	x = o->getx();
					v += o->getv();
					starting = 0;
				 }
			//if(starting==0 && control==1){v+=vrazzi; control=0;}
			if((Taccensioni.size()>0) && (Taccensioni.size()>Nstep) && (T >= Taccensioni[Nstep])) v += vRazzi[Nstep++];
			v += F/M*dT;
			x += v*dT;
		 }
}


void my3main(double a, double b)
{
	double vx=a*10000;
	double vy=b*10000;
	SSolare s;
	crpCeleste sole(1.98e30, V(0., 0.),V(0., 0.));
	crpCeleste mercurio(3.33e23, V(6.97e10, 0.),V(0., 43479.17));
	crpCeleste venere(4.87e24, V(1.09e11, 0.),V(0., 34840.23));
	crpCeleste terra(5.97e24, V(1.52e11, 0.),V(0., 29476.35));
	crpCeleste marte(6.42e23, V(2.49e11, 0.),V(0., 23025.48));
	crpCeleste giove(1.899e27, V(0., 8.16e11),V(-12723.41, 0.));
	crpCeleste saturno(5.685e26, V(1.51e12, 0.),V(0., 9370.69));
	crpCeleste urano(8.68e25, V(3.004e12, 0.),V(0., 6836.23));
	//crpCeleste nettuno(1.02e26,V(4.537e12,0.),V(0.,5478.45));
	Apollo shuttle(110000.);
	shuttle.setv(V(vx,vy));
	shuttle.setOrigin(&terra);
	shuttle.setStartTime(180*86400.);
	shuttle.accendiRazzi(400*86400,V(-5000., 0.));
	s.add(&sole);
	s.add(&mercurio);
	s.add(&venere);
	s.add(&terra);
	s.add(&marte);
	s.add(&giove);
	s.add(&saturno);
	s.add(&urano);
	//s.add(&nettuno);
	s.add(&shuttle);
	//s.setMaxTime(20*3.16e7);
	s.setMaxTime(5*3.16e7);
	s.setDeltaT(86400.);
	c = new TCanvas("c", "Sistema Solare",0, 0, 800, 800);
	c->SetFillColor(1);
	TView *view;
	view = TView::CreateView();
	view->SetRange(-400,-400,-400,400,400,400);
	p = new TPolyMarker*[s.how_many()];
	for (int i=0; i<s.how_many(); i++)
		 {
			if(i==9) ++i;
			p[i] = new TPolyMarker(1);
			p[i]->SetMarkerStyle(20);
			p[i]->SetMarkerColor(10-i);
			p[i]->SetMarkerSize(20);
			p[i]->Draw();
			}
	c->Update();
	s.evolvi(&output);
  c->Update();

}
